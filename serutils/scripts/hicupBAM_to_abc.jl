using ArgParse
using XAM

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--bam"
            help = "HICUP BAM input file"
        "--output", "-o"
            help = "output .tsv file"
        "--resolution", "-r"
            help = "resolution"
            arg_type = Int
            required = true
        "--bin_coord"
            help = "3 columns format: with bin coordinates instead of genomic coordinates (genomic bin instead of chromosome and position)"
            action = :store_true
    end

    return parse_args(s)
end


function check_str2(a)
    return tryparse(Int64, a) !== nothing
end


function extract_bam_5columns(input, output, resolution)
    println("extracting interactions from BAM...")
    reader = open(BAM.Reader, input)

    # get chromosome lengths
    chromosome_lens = Dict(zip(reader.refseqnames, reader.refseqlens))
    # sort chromosome names according to convention (first numbers than letters, weird chromosomes at the end)
    chromosomes = sort(
        reader.refseqnames, 
        by=x -> check_str2(x) ? parse(Int, x) : sum([100 + Int(only(c)) for c in x])
        )

    data = Dict{String, Int64}()
    for record in reader
        chrom1 = BAM.refname(record)
        chrom2 = BAM.nextrefname(record)
        pos1 = Int(BAM.position(record)) ÷ resolution
        pos2 = Int(BAM.nextposition(record)) ÷ resolution
        pos = "$chrom1\t$pos1\t$chrom2\t$pos2"
        if haskey(data, pos)
            data[pos] += 1
        else
            data[pos] = 1
        end
    end
    close(reader)

    println("writing to TSV...")
    open(output, "w") do io
        for n in chromosomes
            l = chromosome_lens[n]
            write(io, "# CHROMOSOME\t$n\t$l\n")
        end
        write(io, "# BIN SIZE\t$resolution\n")
        for (pos, val) in data
            write(io, string(
                "$pos\t$val\n")
                )
        end
    end
end

function extract_bam_3columns(input, output, resolution)
    println("extracting interactions from BAM...")
    reader = open(BAM.Reader, input)

    # get chromosome lengths
    chromosome_lens = Dict(zip(reader.refseqnames, reader.refseqlens))
    # sort chromosome names according to convention (first numbers than letters, weird chromosomes at the end)
    chromosomes = sort(
        reader.refseqnames, 
        by=x -> check_str2(x) ? parse(Int, x) : sum([100 + Int(only(c)) for c in x])
        )
    chromosome_pos = Dict{String, Int64}()
    total = 0
    for chrom in chromosomes
        chromosome_pos[chrom] = total
        total += chromosome_lens[chrom] ÷ resolution + 1
    end

    data = Dict{String, Int32}()
    for record in reader
        chrom1 = BAM.refname(record)
        chrom2 = BAM.nextrefname(record)
        pos1 = Int(BAM.position(record)) ÷ resolution + chromosome_pos[chrom1]
        pos2 = Int(BAM.nextposition(record)) ÷ resolution + chromosome_pos[chrom2]
        pos = "$pos1\t$pos2"
        if haskey(data, pos)
            data[pos] += 1
        else
            data[pos] = 1
        end
    end
    close(reader)

    println("writing to TSV...")
    open(output, "w") do io
        for n in chromosomes
            l = chromosome_lens[n]
            write(io, "# CHROMOSOME\t$n\t$l\n")
        end
        write(io, "# BIN SIZE\t$resolution\n")
        for (pos, val) in data
            write(io, string(
                "$pos\t$val\n")
                )
        end
    end
end

function main()

    args = parse_commandline()

    input = args["bam"]
    output = args["output"]
    resolution = args["resolution"]
    bin_coord = args["bin_coord"]
    if bin_coord
        extract_bam_3columns(input, output, resolution)
    else
        extract_bam_5columns(input, output, resolution)
    end
end

main()