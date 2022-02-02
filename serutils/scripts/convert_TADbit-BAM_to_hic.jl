using ArgParse
using XAM

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--bam"
            help = "TADbit BAM input file"
        "--output", "-o"
            help = "output .hic file"
        "--cpus"
            help = "number of CPUs"
            arg_type = Int
            default = 16
        "--genome"
            help = "Genome name, can be any of hg19 hg38 mm9 mm10"
            required = true
        "--juicer_path"
            help = "path to juicer_tools.jar file"
            default = "/home/fransua/Tools/juicer_tools_1.22.01.jar"
        "--juicer_mem"
            help = "maximum memory for juicer tools"
            default = 30720
    end

    return parse_args(s)
end


function extract_bam(input, tmp)
    reader = open(BAM.Reader, input)
    open(tmp, "w") do io
        for record in reader
            aux = BAM.auxdata(record)
            strand1 = aux["S1"] == 1 ? 0 : 1
            strand2 = aux["S2"] == 1 ? 0 : 1
            strand1 = strand2 = 0
            chrom1 = replace(BAM.refname(record), "MT" => "M")
            chrom2 = replace(BAM.nextrefname(record), "MT" => "M")
            chromnum1 = BAM.refid(record)
            chromnum2 = BAM.nextrefid(record)
            pos1 = Int(BAM.position(record))
            pos2 = Int(BAM.nextposition(record))
            if (chromnum1 == chromnum2 && pos1 > pos2) || chromnum1 > chromnum2
                chrom1, chrom2, strand1, strand2, pos1, pos2 = chrom2, chrom1, strand2, strand1, pos2, pos1
            end
            write(io, string(
                strand1, " ",  # 0 should mean forward
                "chr$chrom1 ",
                pos1, " ", 
                "0 ",  # fragment 1
                strand2, " ",  # 0 should mean forward
                "chr$chrom2 ",
                pos2, " ",
                "1\n")  # fragment 2
                )
        end
        close(reader)
    end
end

function main()

    args = parse_commandline()
    input = args["bam"]
    cpus = args["cpus"]
    output = args["output"]
    genome = args["genome"]
    juicer = args["juicer_path"]
    juicer_mem = args["juicer_mem"]
    tmp = output * "_tmp"
    tmp_sorted = output * "_tmp_sorted"

    println("extracting interactions from BAM...")
    extract_bam(input, tmp)

    println("sorting...")
    run(pipeline(`sort -u -S 5% -k2,2V -k6,6V -k3,3n -k7,7n $tmp`, `uniq`, tmp_sorted))

    println("generating .hic...")
    run(`java -Xms1024m -Xmx$juicer_mem\m -jar $juicer pre -j $cpus $tmp_sorted $output $genome`)

    println("removing temporary...")
    run(`rm -f $tmp $tmp_sorted`)
end

main()