using ArgParse
using FASTX

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--fasta"
            help = "Genome in fasta format"
            required = true
        "--output", "-o"
            help = "output file"
            required = true
        "--renz"
            help = "restriction enzyme name"
            required = true
    end

    return parse_args(s)
end

function read_fasta(fasta)
    genome = Dict{String, Array}()
    open(FASTA.Reader, fasta) do fasta_handler
        for record in fasta_handler
            genome[FASTA.identifier(record)] = [v[1] for v in findall(r"GATC", uppercase(FASTA.sequence(String, record)))]
            genome[FASTA.identifier(record)] *= []
        end
    end
    return genome
end


function main()
    args = parse_commandline()
    renz = args["renz"]
    fasta = args["fasta"]
    output = args["output"]

    genome = read_fasta(fasta)

    open(output, "w") do out_handler
        for (chrom, sequence) in genome
            write(out_handler, chrom * " " * (join(sequence, " ")) * "\n")
        end
    end

end

RESTRICTION_ENZYMES = Dict{String, String}(
    "aani" => "TTA|TAA"                     ,
    "aari" => "CACCTGC|"                    ,
    "aasi" => "GACNNNN|NNGTC"               ,
    "aatii" => "GACGT|C"                     ,
    "abasi" => "C|"                          ,
    "absi" => "CC|TCGAGG"                   ,
    "acc16i" => "TGC|GCA"                     ,
    "acc36i" => "ACCTGC|"                     ,
    "acc65i" => "G|GTACC"                     ,
    "accb1i" => "G|GYRCC"                     ,
    "accb7i" => "CCANNNN|NTGG"                ,
    "accbsi" => "CCG|CTC"                     ,
    "acci" => "GT|MKAC"                     ,
    "accii" => "CG|CG"                       ,
    "acciii" => "T|CCGGA"                     ,
    "aceiii" => "CAGCTC|"                     ,
    "acii" => "C|CGC"                       ,
    "acli" => "AA|CGTT"                     ,
    "aclwi" => "GGATC|"                      ,
    "acoi" => "Y|GGCCR"                     ,
    "acsi" => "R|AATTY"                     ,
    "acui" => "CTGAAG|"                     ,
    "acvi" => "CAC|GTG"                     ,
    "acyi" => "GR|CGYC"                     ,
    "adei" => "CACNNN|GTG"                  ,
    "afai" => "GT|AC"                       ,
    "afei" => "AGC|GCT"                     ,
    "afii" => "CCNNNNN|NNGG"                ,
    "aflii" => "C|TTAAG"                     ,
    "afliii" => "A|CRYGT"                     ,
    "agei" => "A|CCGGT"                     ,
    "agsi" => "TTS|AA"                      ,
    "ahaiii" => "TTT|AAA"                     ,
    "ahdi" => "GACNNN|NNGTC"                ,
    "ahli" => "A|CTAGT"                     ,
    "ajii" => "CAC|GTC"                     ,
    "ajni" => "|CCWGG"                      ,
    "ajui" => "GAANNNN|NNNTTGG"             ,
    "alei" => "CACNN|NNGTG"                 ,
    "alfi" => "GC|ANNNNNNTGC"               ,
    "aloi" => "GAACNN|NNNNTCC"              ,
    "alubi" => "AG|CT"                       ,
    "alui" => "AG|CT"                       ,
    "alw21i" => "GWGCW|C"                     ,
    "alw26i" => "GTCTC|"                      ,
    "alw44i" => "G|TGCAC"                     ,
    "alwfi" => "GAAAYNNNNNRTG|GAAAYNNNNNRTG" ,
    "alwi" => "GGATC|"                      ,
    "alwni" => "CAGNNN|CTG"                  ,
    "ama87i" => "C|YCGRG"                     ,
    "aor13hi" => "T|CCGGA"                     ,
    "aor51hi" => "AGC|GCT"                     ,
    "aoxi" => "|GGCC"                       ,
    "apabi" => "GCANNNNN|TGC"                ,
    "apai" => "GGGCC|C"                     ,
    "apali" => "G|TGCAC"                     ,
    "apeki" => "G|CWGC"                      ,
    "apoi" => "R|AATTY"                     ,
    "apypi" => "ATCGAC|"                     ,
    "aquii" => "GCCGNAC|"                    ,
    "aquiii" => "GAGGAG|"                     ,
    "aquiv" => "GRGGAAG|"                    ,
    "arsi" => "GACNN|NNNNTTYG"              ,
    "asci" => "GG|CGCGCC"                   ,
    "asei" => "AT|TAAT"                     ,
    "asi256i" => "G|ATC"                       ,
    "asigi" => "A|CCGGT"                     ,
    "asisi" => "GCGAT|CGC"                   ,
    "asp700i" => "GAANN|NNTTC"                 ,
    "asp718i" => "G|GTACC"                     ,
    "aspa2i" => "C|CTAGG"                     ,
    "aspbhi" => "YSCNS|"                      ,
    "asplei" => "GCG|C"                       ,
    "asps9i" => "G|GNCC"                      ,
    "assi" => "AGT|ACT"                     ,
    "asuc2i" => "CC|SGG"                      ,
    "asuhpi" => "GGTGA|"                      ,
    "asui" => "G|GNCC"                      ,
    "asuii" => "TT|CGAA"                     ,
    "asunhi" => "G|CTAGC"                     ,
    "avai" => "C|YCGRG"                     ,
    "avaii" => "G|GWCC"                      ,
    "avaiii" => "ATGCAT|ATGCAT"               ,
    "avrii" => "C|CTAGG"                     ,
    "axyi" => "CC|TNAGG"                    ,
    "baegi" => "GKGCM|C"                     ,
    "baei" => "A|CNNNNGTAYC"                ,
    "bali" => "TGG|CCA"                     ,
    "bamhi" => "G|GATCC"                     ,
    "bani" => "G|GYRCC"                     ,
    "banii" => "GRGCY|C"                     ,
    "bari" => "GAAGNN|NNNNTAC"              ,
    "basi" => "CCANNNN|NTGG"                ,
    "baui" => "C|ACGAG"                     ,
    "bbr7i" => "GAAGAC|"                     ,
    "bbrpi" => "CAC|GTG"                     ,
    "bbsi" => "GAAGAC|"                     ,
    "bbv12i" => "GWGCW|C"                     ,
    "bbvci" => "CC|TCAGC"                    ,
    "bbvi" => "GCAGC|"                      ,
    "bbvii" => "GAAGAC|"                     ,
    "bcci" => "CCATC|"                      ,
    "bce83i" => "CTTGAG|"                     ,
    "bceai" => "ACGGC|"                      ,
    "bcefi" => "ACGGC|"                      ,
    "bcgi" => "CG|ANNNNNNTGC"               ,
    "bcit130i" => "CC|WGG"                      ,
    "bcivi" => "GTATCC|"                     ,
    "bcli" => "T|GATCA"                     ,
    "bcni" => "CC|SGG"                      ,
    "bcodi" => "GTCTC|"                      ,
    "bcui" => "A|CTAGT"                     ,
    "bdai" => "TG|ANNNNNNTCA"               ,
    "beti" => "W|CCGGW"                     ,
    "bfai" => "C|TAG"                       ,
    "bfii" => "ACTGGG|"                     ,
    "bfmi" => "C|TRYAG"                     ,
    "bfoi" => "RGCGC|Y"                     ,
    "bfri" => "C|TTAAG"                     ,
    "bfuai" => "ACCTGC|"                     ,
    "bfuci" => "|GATC"                       ,
    "bfui" => "GTATCC|"                     ,
    "bgli" => "GCCNNNN|NGGC"                ,
    "bglii" => "A|GATCT"                     ,
    "bini" => "GGATC|"                      ,
    "bisi" => "GC|NGC"                      ,
    "blni" => "C|CTAGG"                     ,
    "blpi" => "GC|TNAGC"                    ,
    "blsi" => "GCN|GC"                      ,
    "bmcai" => "AGT|ACT"                     ,
    "bme1390i" => "CC|NGG"                      ,
    "bme18i" => "G|GWCC"                      ,
    "bmedi" => "C|"                          ,
    "bmeri" => "GACNNN|NNGTC"                ,
    "bmet110i" => "C|YCGRG"                     ,
    "bmgbi" => "CAC|GTC"                     ,
    "bmgi" => "GKGCCC|GKGCCC"               ,
    "bmgt120i" => "G|GNCC"                      ,
    "bmii" => "GGN|NCC"                     ,
    "bmrfi" => "CC|NGG"                      ,
    "bmri" => "ACTGGG|"                     ,
    "bmsi" => "GCATC|"                      ,
    "bmti" => "GCTAG|C"                     ,
    "bmui" => "ACTGGG|"                     ,
    "boxi" => "GACNN|NNGTC"                 ,
    "bpii" => "GAAGAC|"                     ,
    "bpli" => "GAG|NNNNNCTC"                ,
    "bpmi" => "CTGGAG|"                     ,
    "bpu10i" => "CC|TNAGC"                    ,
    "bpu1102i" => "GC|TNAGC"                    ,
    "bpu14i" => "TT|CGAA"                     ,
    "bpuei" => "CTTGAG|"                     ,
    "bpumi" => "CC|SGG"                      ,
    "bpvui" => "CGAT|CG"                     ,
    "bsa29i" => "AT|CGAT"                     ,
    "bsaai" => "YAC|GTR"                     ,
    "bsabi" => "GATNN|NNATC"                 ,
    "bsahi" => "GR|CGYC"                     ,
    "bsai" => "GGTCTC|"                     ,
    "bsaji" => "C|CNNGG"                     ,
    "bsawi" => "W|CCGGW"                     ,
    "bsaxi" => "AC|NNNNNCTCC"                ,
    "bsbi" => "CAACAC|"                     ,
    "bsc4i" => "CCNNNNN|NNGG"                ,
    "bscai" => "GCATC|"                      ,
    "bscgi" => "CCCGT|CCCGT"                 ,
    "bse118i" => "R|CCGGY"                     ,
    "bse1i" => "ACTGG|"                      ,
    "bse21i" => "CC|TNAGG"                    ,
    "bse3di" => "GCAATG|"                     ,
    "bse8i" => "GATNN|NNATC"                 ,
    "bseai" => "T|CCGGA"                     ,
    "bsebi" => "CC|WGG"                      ,
    "bseci" => "AT|CGAT"                     ,
    "bsedi" => "C|CNNGG"                     ,
    "bsegi" => "GGATG|"                      ,
    "bseji" => "GATNN|NNATC"                 ,
    "bseli" => "CCNNNNN|NNGG"                ,
    "bsemi" => "GCAATG|"                     ,
    "bsemii" => "CTCAG|"                      ,
    "bseni" => "ACTGG|"                      ,
    "bsepi" => "G|CGCGC"                     ,
    "bseri" => "GAGGAG|"                     ,
    "bsesi" => "GKGCM|C"                     ,
    "bsex3i" => "C|GGCCG"                     ,
    "bsexi" => "GCAGC|"                      ,
    "bseyi" => "C|CCAGC"                     ,
    "bsgi" => "GTGCAG|"                     ,
    "bsh1236i" => "CG|CG"                       ,
    "bsh1285i" => "CGRY|CG"                     ,
    "bshfi" => "GG|CC"                       ,
    "bshni" => "G|GYRCC"                     ,
    "bshti" => "A|CCGGT"                     ,
    "bshvi" => "AT|CGAT"                     ,
    "bsiei" => "CGRY|CG"                     ,
    "bsihkai" => "GWGCW|C"                     ,
    "bsihkci" => "C|YCGRG"                     ,
    "bsii" => "C|ACGAG"                     ,
    "bsisi" => "C|CGG"                       ,
    "bsiwi" => "C|GTACG"                     ,
    "bsiyi" => "CCNNNNN|NNGG"                ,
    "bslfi" => "GGGAC|"                      ,
    "bsli" => "CCNNNNN|NNGG"                ,
    "bsmai" => "GTCTC|"                      ,
    "bsmbi" => "CGTCTC|"                     ,
    "bsmfi" => "GGGAC|"                      ,
    "bsmi" => "GAATGC|"                     ,
    "bsni" => "GG|CC"                       ,
    "bso31i" => "GGTCTC|"                     ,
    "bsobi" => "C|YCGRG"                     ,
    "bsp119i" => "TT|CGAA"                     ,
    "bsp120i" => "G|GGCCC"                     ,
    "bsp1286i" => "GDGCH|C"                     ,
    "bsp13i" => "T|CCGGA"                     ,
    "bsp1407i" => "T|GTACA"                     ,
    "bsp143i" => "|GATC"                       ,
    "bsp1720i" => "GC|TNAGC"                    ,
    "bsp19i" => "C|CATGG"                     ,
    "bsp24i" => "GACN|NNNNNTGG"               ,
    "bsp68i" => "TCG|CGA"                     ,
    "bspaci" => "C|CGC"                       ,
    "bspcni" => "CTCAG|"                      ,
    "bspd6i" => "GACTC|"                      ,
    "bspdi" => "AT|CGAT"                     ,
    "bspei" => "T|CCGGA"                     ,
    "bspfni" => "CG|CG"                       ,
    "bspgi" => "CTGGAC|CTGGAC"               ,
    "bsphi" => "T|CATGA"                     ,
    "bspli" => "GGN|NCC"                     ,
    "bsplu11i" => "A|CATGT"                     ,
    "bspmi" => "ACCTGC|"                     ,
    "bspmii" => "T|CCGGA"                     ,
    "bspnci" => "CCAGA|CCAGA"                 ,
    "bspoi" => "GCTAG|C"                     ,
    "bsppi" => "GGATC|"                      ,
    "bspqi" => "GCTCTTC|"                    ,
    "bspt104i" => "TT|CGAA"                     ,
    "bspt107i" => "G|GYRCC"                     ,
    "bspti" => "C|TTAAG"                     ,
    "bsrbi" => "CCG|CTC"                     ,
    "bsrdi" => "GCAATG|"                     ,
    "bsrfi" => "R|CCGGY"                     ,
    "bsrgi" => "T|GTACA"                     ,
    "bsri" => "ACTGG|"                      ,
    "bsrsi" => "ACTGG|"                      ,
    "bssai" => "R|CCGGY"                     ,
    "bsseci" => "C|CNNGG"                     ,
    "bsshii" => "G|CGCGC"                     ,
    "bsski" => "|CCNGG"                      ,
    "bssmi" => "|GATC"                       ,
    "bssnai" => "GTA|TAC"                     ,
    "bssni" => "GR|CGYC"                     ,
    "bsssi" => "C|ACGAG"                     ,
    "bsst1i" => "C|CWWGG"                     ,
    "bst1107i" => "GTA|TAC"                     ,
    "bst2bi" => "C|ACGAG"                     ,
    "bst2ui" => "CC|WGG"                      ,
    "bst4ci" => "ACN|GT"                      ,
    "bst6i" => "CTCTTC|"                     ,
    "bstaci" => "GR|CGYC"                     ,
    "bstafi" => "C|TTAAG"                     ,
    "bstapi" => "GCANNNN|NTGC"                ,
    "bstaui" => "T|GTACA"                     ,
    "bstbai" => "YAC|GTR"                     ,
    "bstbi" => "TT|CGAA"                     ,
    "bstc8i" => "GCN|NGC"                     ,
    "bstdei" => "C|TNAG"                      ,
    "bstdsi" => "C|CRYGG"                     ,
    "bsteii" => "G|GTNACC"                    ,
    "bsteni" => "CCTNN|NNNAGG"                ,
    "bstf5i" => "GGATG|"                      ,
    "bstfni" => "CG|CG"                       ,
    "bsth2i" => "RGCGC|Y"                     ,
    "bsthhi" => "GCG|C"                       ,
    "bstkti" => "GAT|C"                       ,
    "bstmai" => "GTCTC|"                      ,
    "bstmbi" => "|GATC"                       ,
    "bstmci" => "CGRY|CG"                     ,
    "bstmwi" => "GCNNNNN|NNGC"                ,
    "bstni" => "CC|WGG"                      ,
    "bstnsi" => "RCATG|Y"                     ,
    "bstoi" => "CC|WGG"                      ,
    "bstpai" => "GACNN|NNGTC"                 ,
    "bstpi" => "G|GTNACC"                    ,
    "bstsci" => "|CCNGG"                      ,
    "bstsfi" => "C|TRYAG"                     ,
    "bstsli" => "GKGCM|C"                     ,
    "bstsni" => "TAC|GTA"                     ,
    "bstui" => "CG|CG"                       ,
    "bstv1i" => "GCAGC|"                      ,
    "bstv2i" => "GAAGAC|"                     ,
    "bstx2i" => "R|GATCY"                     ,
    "bstxi" => "CCANNNNN|NTGG"               ,
    "bstyi" => "R|GATCY"                     ,
    "bstz17i" => "GTA|TAC"                     ,
    "bstzi" => "C|GGCCG"                     ,
    "bsu15i" => "AT|CGAT"                     ,
    "bsu36i" => "CC|TNAGG"                    ,
    "bsui" => "GTATCC|"                     ,
    "bsuri" => "GG|CC"                       ,
    "btgi" => "C|CRYGG"                     ,
    "btgzi" => "GCGATG|"                     ,
    "bthci" => "GCNG|C"                      ,
    "btri" => "CAC|GTC"                     ,
    "btsci" => "GGATG|"                      ,
    "btsi" => "GCAGTG|"                     ,
    "btsimuti" => "CAGTG|"                      ,
    "btumi" => "TCG|CGA"                     ,
    "bvei" => "ACCTGC|"                     ,
    "cac8i" => "GCN|NGC"                     ,
    "caii" => "CAGNNN|CTG"                  ,
    "cauii" => "CC|SGG"                      ,
    "cchii" => "GGARGA|"                     ,
    "cchiii" => "CCCAAG|"                     ,
    "ccii" => "T|CATGA"                     ,
    "ccini" => "GC|GGCCGC"                   ,
    "cdi630v" => "CAAAAA|CAAAAA"               ,
    "cdii" => "CATC|G"                      ,
    "cdpi" => "GCGGAG|"                     ,
    "cfoi" => "GCG|C"                       ,
    "cfr10i" => "R|CCGGY"                     ,
    "cfr13i" => "G|GNCC"                      ,
    "cfr42i" => "CCGC|GG"                     ,
    "cfr9i" => "C|CCGGG"                     ,
    "cfri" => "Y|GGCCR"                     ,
    "cgl13032i" => "GGCGCA|GGCGCA"               ,
    "cgl13032ii" => "ACGABGG|ACGABGG"             ,
    "chai" => "GATC|"                       ,
    "cjefiii" => "GCAAGG|GCAAGG"               ,
    "cjefv" => "GGRCA|GGRCA"                 ,
    "cjei" => "CCA|NNNNNNGT"                ,
    "cjenii" => "GAGNNNNNGT|GAGNNNNNGT"       ,
    "cjeniii" => "GKAAYG|"                     ,
    "cjep659iv" => "CACNNNNNNNGAA|CACNNNNNNNGAA" ,
    "cjepi" => "CCANN|NNNNNTC"               ,
    "cjui" => "CAYNNNNNRTG|CAYNNNNNRTG"     ,
    "cjuii" => "CAYNNNNNCTC|CAYNNNNNCTC"     ,
    "clai" => "AT|CGAT"                     ,
    "cpoi" => "CG|GWCCG"                    ,
    "csei" => "GACGC|"                      ,
    "csii" => "A|CCWGGT"                    ,
    "csp6i" => "G|TAC"                       ,
    "cspai" => "A|CCGGT"                     ,
    "cspci" => "C|AANNNNNGTGG"               ,
    "cspi" => "CG|GWCCG"                    ,
    "cstmi" => "AAGGAG|"                     ,
    "cviaii" => "C|ATG"                       ,
    "cviji" => "RG|CY"                       ,
    "cviki_1" => "RG|CY"                       ,
    "cviqi" => "G|TAC"                       ,
    "cviri" => "TG|CA"                       ,
    "ddei" => "C|TNAG"                      ,
    "dini" => "GGC|GCC"                     ,
    "dpni" => "GA|TC"                       ,
    "dpnii" => "|GATC"                       ,
    "drai" => "TTT|AAA"                     ,
    "draii" => "RG|GNCCY"                    ,
    "draiii" => "CACNNN|GTG"                  ,
    "drari" => "CAAGNAC|"                    ,
    "drdi" => "GACNNNN|NNGTC"               ,
    "drdii" => "GAACCA|GAACCA"               ,
    "drii" => "GACNNN|NNGTC"                ,
    "dsai" => "C|CRYGG"                     ,
    "dsedi" => "GACNNNN|NNGTC"               ,
    "eaei" => "Y|GGCCR"                     ,
    "eagi" => "C|GGCCG"                     ,
    "eam1104i" => "CTCTTC|"                     ,
    "eam1105i" => "GACNNN|NNGTC"                ,
    "eari" => "CTCTTC|"                     ,
    "ecii" => "GGCGGA|"                     ,
    "ecl136ii" => "GAG|CTC"                     ,
    "eclxi" => "C|GGCCG"                     ,
    "eco105i" => "TAC|GTA"                     ,
    "eco130i" => "C|CWWGG"                     ,
    "eco147i" => "AGG|CCT"                     ,
    "eco24i" => "GRGCY|C"                     ,
    "eco31i" => "GGTCTC|"                     ,
    "eco32i" => "GAT|ATC"                     ,
    "eco47i" => "G|GWCC"                      ,
    "eco47iii" => "AGC|GCT"                     ,
    "eco52i" => "C|GGCCG"                     ,
    "eco53ki" => "GAG|CTC"                     ,
    "eco57i" => "CTGAAG|"                     ,
    "eco57mi" => "CTGRAG|"                     ,
    "eco72i" => "CAC|GTG"                     ,
    "eco81i" => "CC|TNAGG"                    ,
    "eco88i" => "C|YCGRG"                     ,
    "eco91i" => "G|GTNACC"                    ,
    "ecohi" => "|CCSGG"                      ,
    "ecoicri" => "GAG|CTC"                     ,
    "econi" => "CCTNN|NNNAGG"                ,
    "ecoo109i" => "RG|GNCCY"                    ,
    "ecoo65i" => "G|GTNACC"                    ,
    "ecori" => "G|AATTC"                     ,
    "ecorii" => "|CCWGG"                      ,
    "ecorv" => "GAT|ATC"                     ,
    "ecot14i" => "C|CWWGG"                     ,
    "ecot22i" => "ATGCA|T"                     ,
    "ecot38i" => "GRGCY|C"                     ,
    "egei" => "GGC|GCC"                     ,
    "ehei" => "GGC|GCC"                     ,
    "erhi" => "C|CWWGG"                     ,
    "esabc3i" => "TC|GA"                       ,
    "esassi" => "GACCAC|GACCAC"               ,
    "esp3i" => "CGTCTC|"                     ,
    "espi" => "GC|TNAGC"                    ,
    "faei" => "CATG|"                       ,
    "faii" => "YA|TR"                       ,
    "fali" => "AAG|NNNNNCTT"                ,
    "faqi" => "GGGAC|"                      ,
    "fati" => "|CATG"                       ,
    "faui" => "CCCGC|"                      ,
    "faundi" => "CA|TATG"                     ,
    "fbai" => "T|GATCA"                     ,
    "fbli" => "GT|MKAC"                     ,
    "fini" => "GGGAC|GGGAC"                 ,
    "fmui" => "GGNC|C"                      ,
    "fnu4hi" => "GC|NGC"                      ,
    "fnudii" => "CG|CG"                       ,
    "foki" => "GGATG|"                      ,
    "frioi" => "GRGCY|C"                     ,
    "fsei" => "GGCCGG|CC"                   ,
    "fsp4hi" => "GC|NGC"                      ,
    "fspai" => "RTGC|GCAY"                   ,
    "fspbi" => "C|TAG"                       ,
    "fspei" => "CC|"                         ,
    "fspi" => "TGC|GCA"                     ,
    "gaut27i" => "CGCGCAGG|CGCGCAGG"           ,
    "gdiii" => "C|GGCCR"                     ,
    "glai" => "GC|GC"                       ,
    "glui" => "GC|NGC"                      ,
    "gsai" => "CCCAG|C"                     ,
    "gsui" => "CTGGAG|"                     ,
    "haei" => "WGG|CCW"                     ,
    "haeii" => "RGCGC|Y"                     ,
    "haeiii" => "GG|CC"                       ,
    "hapii" => "C|CGG"                       ,
    "hauii" => "TGGCCA|"                     ,
    "hgai" => "GACGC|"                      ,
    "hgiai" => "GWGCW|C"                     ,
    "hgici" => "G|GYRCC"                     ,
    "hgieii" => "ACCNNNNNNGGT|ACCNNNNNNGGT"   ,
    "hgijii" => "GRGCY|C"                     ,
    "hhai" => "GCG|C"                       ,
    "hin1i" => "GR|CGYC"                     ,
    "hin1ii" => "CATG|"                       ,
    "hin4i" => "GAY|NNNNNVTC"                ,
    "hin4ii" => "CCTTC|"                      ,
    "hin6i" => "G|CGC"                       ,
    "hinp1i" => "G|CGC"                       ,
    "hincii" => "GTY|RAC"                     ,
    "hindii" => "GTY|RAC"                     ,
    "hindiii" => "A|AGCTT"                     ,
    "hinfi" => "G|ANTC"                      ,
    "hpai" => "GTT|AAC"                     ,
    "hpaii" => "C|CGG"                       ,
    "hphi" => "GGTGA|"                      ,
    "hpy166ii" => "GTN|NAC"                     ,
    "hpy178iii" => "TC|NNGA"                     ,
    "hpy188i" => "TCN|GA"                      ,
    "hpy188iii" => "TC|NNGA"                     ,
    "hpy8i" => "GTN|NAC"                     ,
    "hpy99i" => "CGWCG|"                      ,
    "hpy99xiii" => "GCCTA|GCCTA"                 ,
    "hpy99xiv" => "GGWTAA|GGWTAA"               ,
    "hpyav" => "CCTTC|"                      ,
    "hpych4iii" => "ACN|GT"                      ,
    "hpych4iv" => "A|CGT"                       ,
    "hpych4v" => "TG|CA"                       ,
    "hpyf10vi" => "GCNNNNN|NNGC"                ,
    "hpyf3i" => "C|TNAG"                      ,
    "hpyse526i" => "A|CGT"                       ,
    "hsp92i" => "GR|CGYC"                     ,
    "hsp92ii" => "CATG|"                       ,
    "hspai" => "G|CGC"                       ,
    "jma19592i" => "GTATNAC|GTATNAC"             ,
    "kasi" => "G|GCGCC"                     ,
    "kfli" => "GG|GWCCC"                    ,
    "kpn2i" => "T|CCGGA"                     ,
    "kpni" => "GGTAC|C"                     ,
    "kroi" => "G|CCGGC"                     ,
    "ksp22i" => "T|GATCA"                     ,
    "ksp632i" => "CTCTTC|"                     ,
    "kspai" => "GTT|AAC"                     ,
    "kspi" => "CCGC|GG"                     ,
    "kzo9i" => "|GATC"                       ,
    "lgui" => "GCTCTTC|"                    ,
    "lpni" => "RGC|GCY"                     ,
    "lpnpi" => "CCDG|"                       ,
    "lsp1109i" => "GCAGC|"                      ,
    "lwei" => "GCATC|"                      ,
    "mabi" => "A|CCWGGT"                    ,
    "maei" => "C|TAG"                       ,
    "maeii" => "A|CGT"                       ,
    "maeiii" => "|GTNAC"                      ,
    "mali" => "GA|TC"                       ,
    "maqi" => "CRTTGAC|"                    ,
    "maubi" => "CG|CGCGCG"                   ,
    "mbii" => "CCG|CTC"                     ,
    "mboi" => "|GATC"                       ,
    "mboii" => "GAAGA|"                      ,
    "mcati" => "GCGC|GC"                     ,
    "mcri" => "CGRY|CG"                     ,
    "mfei" => "C|AATTG"                     ,
    "mfli" => "R|GATCY"                     ,
    "mhli" => "GDGCH|C"                     ,
    "mjaiv" => "GTNNAC|GTNNAC"               ,
    "mkadii" => "GAGAYGT|GAGAYGT"             ,
    "mlsi" => "TGG|CCA"                     ,
    "mluci" => "|AATT"                       ,
    "mlui" => "A|CGCGT"                     ,
    "mluni" => "TGG|CCA"                     ,
    "mly113i" => "GG|CGCC"                     ,
    "mlyi" => "GAGTC|"                      ,
    "mmei" => "TCCRAC|"                     ,
    "mnli" => "CCTC|"                       ,
    "mph1103i" => "ATGCA|T"                     ,
    "mrei" => "CG|CCGGCG"                   ,
    "mroi" => "T|CCGGA"                     ,
    "mroni" => "G|CCGGC"                     ,
    "mroxi" => "GAANN|NNTTC"                 ,
    "msci" => "TGG|CCA"                     ,
    "msei" => "T|TAA"                       ,
    "msli" => "CAYNN|NNRTG"                 ,
    "msp20i" => "TGG|CCA"                     ,
    "mspa1i" => "CMG|CKG"                     ,
    "mspci" => "C|TTAAG"                     ,
    "mspi" => "C|CGG"                       ,
    "mspji" => "CNNR|"                       ,
    "mspr9i" => "CC|NGG"                      ,
    "mssi" => "GTTT|AAAC"                   ,
    "msti" => "TGC|GCA"                     ,
    "muni" => "C|AATTG"                     ,
    "mva1269i" => "GAATGC|"                     ,
    "mvai" => "CC|WGG"                      ,
    "mvni" => "CG|CG"                       ,
    "mvri" => "CGAT|CG"                     ,
    "mwoi" => "GCNNNNN|NNGC"                ,
    "naei" => "GCC|GGC"                     ,
    "nari" => "GG|CGCC"                     ,
    "ncii" => "CC|SGG"                      ,
    "ncoi" => "C|CATGG"                     ,
    "ndei" => "CA|TATG"                     ,
    "ndeii" => "|GATC"                       ,
    "ngoaviii" => "|GACNNNNNTGA"                ,
    "ngomiv" => "G|CCGGC"                     ,
    "nhaxi" => "CAAGRAG|CAAGRAG"             ,
    "nhei" => "G|CTAGC"                     ,
    "nlaci" => "CATCAC|"                     ,
    "nlaiii" => "CATG|"                       ,
    "nlaiv" => "GGN|NCC"                     ,
    "nli3877i" => "CYCGR|G"                     ,
    "nmeaiii" => "GCCGAG|"                     ,
    "nmedi" => "|RCCGGY"                     ,
    "nmuci" => "|GTSAC"                      ,
    "noti" => "GC|GGCCGC"                   ,
    "nrui" => "TCG|CGA"                     ,
    "nsbi" => "TGC|GCA"                     ,
    "nsii" => "ATGCA|T"                     ,
    "nspbii" => "CMG|CKG"                     ,
    "nspi" => "RCATG|Y"                     ,
    "nspv" => "TT|CGAA"                     ,
    "olii" => "CACNN|NNGTG"                 ,
    "pabi" => "GTA|C"                       ,
    "paci" => "TTAAT|TAA"                   ,
    "paei" => "GCATG|C"                     ,
    "paer7i" => "C|TCGAG"                     ,
    "pagi" => "T|CATGA"                     ,
    "palai" => "GG|CGCGCC"                   ,
    "pasi" => "CC|CWGGG"                    ,
    "paui" => "G|CGCGC"                     ,
    "pcei" => "AGG|CCT"                     ,
    "pcii" => "A|CATGT"                     ,
    "pcisi" => "GCTCTTC|"                    ,
    "pcsi" => "WCGNNNN|NNNCGW"              ,
    "pcti" => "GAATGC|"                     ,
    "pdii" => "GCC|GGC"                     ,
    "pdmi" => "GAANN|NNTTC"                 ,
    "peni" => "GCAGT|GCAGT"                 ,
    "pfei" => "G|AWTC"                      ,
    "pfl1108i" => "TCGTAG|TCGTAG"               ,
    "pfl23ii" => "C|GTACG"                     ,
    "pflfi" => "GACN|NNGTC"                  ,
    "pflmi" => "CCANNNN|NTGG"                ,
    "pfoi" => "T|CCNGGA"                    ,
    "pinai" => "A|CCGGT"                     ,
    "pladi" => "CATCAG|"                     ,
    "ple19i" => "CGAT|CG"                     ,
    "plei" => "GAGTC|"                      ,
    "pluti" => "GGCGC|C"                     ,
    "pmaci" => "CAC|GTG"                     ,
    "pmei" => "GTTT|AAAC"                   ,
    "pmli" => "CAC|GTG"                     ,
    "ppii" => "GAACN|NNNNCTC"               ,
    "ppsi" => "GAGTC|"                      ,
    "ppu10i" => "A|TGCAT"                     ,
    "ppu21i" => "YAC|GTR"                     ,
    "ppumi" => "RG|GWCCY"                    ,
    "psci" => "A|CATGT"                     ,
    "pshai" => "GACNN|NNGTC"                 ,
    "pshbi" => "AT|TAAT"                     ,
    "psii" => "TTA|TAA"                     ,
    "psp03i" => "GGWC|C"                      ,
    "psp124bi" => "GAGCT|C"                     ,
    "psp1406i" => "AA|CGTT"                     ,
    "psp5ii" => "RG|GWCCY"                    ,
    "psp6i" => "|CCWGG"                      ,
    "pspci" => "CAC|GTG"                     ,
    "pspei" => "G|GTNACC"                    ,
    "pspgi" => "|CCWGG"                      ,
    "pspli" => "C|GTACG"                     ,
    "pspn4i" => "GGN|NCC"                     ,
    "pspomi" => "G|GGCCC"                     ,
    "pspomii" => "CGCCCAR|"                    ,
    "psppi" => "G|GNCC"                      ,
    "pspppi" => "RG|GWCCY"                    ,
    "psppri" => "CCYCAG|"                     ,
    "pspxi" => "VC|TCGAGB"                   ,
    "psri" => "GAACNN|NNNNTAC"              ,
    "pssi" => "RGGNC|CY"                    ,
    "psti" => "CTGCA|G"                     ,
    "pstni" => "CAGNNN|CTG"                  ,
    "psui" => "R|GATCY"                     ,
    "psyi" => "GACN|NNGTC"                  ,
    "ptei" => "G|CGCGC"                     ,
    "pvui" => "CGAT|CG"                     ,
    "pvuii" => "CAG|CTG"                     ,
    "r2_bcesiv" => "|GCAGC"                      ,
    "rcei" => "CATCGAC|"                    ,
    "rdegbi" => "CCGCAG|CCGCAG"               ,
    "rdegbii" => "ACCCAG|"                     ,
    "rdegbiii" => "|TGRYCA"                     ,
    "rflfiii" => "CGCCAG|CGCCAG"               ,
    "rgai" => "GCGAT|CGC"                   ,
    "rigi" => "GGCCGG|CC"                   ,
    "rlai" => "VCW|VCW"                     ,
    "rleai" => "CCCACA|"                     ,
    "rpab5i" => "CGRGGAC|"                    ,
    "rpabi" => "CCCGCAG|"                    ,
    "rpai" => "GTYGGAG|"                    ,
    "rpati" => "GRTGGAG|GRTGGAG"             ,
    "rrui" => "TCG|CGA"                     ,
    "rsai" => "GT|AC"                       ,
    "rsani" => "G|TAC"                       ,
    "rsei" => "CAYNN|NNRTG"                 ,
    "rsr2i" => "CG|GWCCG"                    ,
    "rsrii" => "CG|GWCCG"                    ,
    "saci" => "GAGCT|C"                     ,
    "sacii" => "CCGC|GG"                     ,
    "sali" => "G|TCGAC"                     ,
    "sandi" => "GG|GWCCC"                    ,
    "sapi" => "GCTCTTC|"                    ,
    "saqai" => "T|TAA"                       ,
    "sati" => "GC|NGC"                      ,
    "sau3ai" => "|GATC"                       ,
    "sau96i" => "G|GNCC"                      ,
    "saui" => "CC|TNAGG"                    ,
    "sbfi" => "CCTGCA|GG"                   ,
    "scai" => "AGT|ACT"                     ,
    "schi" => "GAGTC|"                      ,
    "scii" => "CTC|GAG"                     ,
    "scrfi" => "CC|NGG"                      ,
    "sdai" => "CCTGCA|GG"                   ,
    "sdeai" => "CAGRAG|"                     ,
    "sdeosi" => "|GACNNNNRTGA"                ,
    "sdui" => "GDGCH|C"                     ,
    "seci" => "C|CNNGG"                     ,
    "seli" => "|CGCG"                       ,
    "seti" => "ASST|"                       ,
    "sexai" => "A|CCWGGT"                    ,
    "sfaai" => "GCGAT|CGC"                   ,
    "sfani" => "GCATC|"                      ,
    "sfci" => "C|TRYAG"                     ,
    "sfei" => "C|TRYAG"                     ,
    "sfii" => "GGCCNNNN|NGGCC"              ,
    "sfoi" => "GGC|GCC"                     ,
    "sfr274i" => "C|TCGAG"                     ,
    "sfr303i" => "CCGC|GG"                     ,
    "sfui" => "TT|CGAA"                     ,
    "sgei" => "CNNG|"                       ,
    "sgfi" => "GCGAT|CGC"                   ,
    "sgrai" => "CR|CCGGYG"                   ,
    "sgrbi" => "CCGC|GG"                     ,
    "sgrdi" => "CG|TCGACG"                   ,
    "sgrti" => "CCDS|"                       ,
    "sgsi" => "GG|CGCGCC"                   ,
    "simi" => "GG|GTC"                      ,
    "slai" => "C|TCGAG"                     ,
    "smai" => "CCC|GGG"                     ,
    "smii" => "ATTT|AAAT"                   ,
    "smimi" => "CAYNN|NNRTG"                 ,
    "smli" => "C|TYRAG"                     ,
    "smoi" => "C|TYRAG"                     ,
    "snabi" => "TAC|GTA"                     ,
    "snai" => "GTATAC|GTATAC"               ,
    "sno506i" => "GGCCGAG|GGCCGAG"             ,
    "spei" => "A|CTAGT"                     ,
    "sphi" => "GCATG|C"                     ,
    "spli" => "C|GTACG"                     ,
    "spodi" => "GCGGRAG|GCGGRAG"             ,
    "srfi" => "GCCC|GGGC"                   ,
    "sse232i" => "CG|CCGGCG"                   ,
    "sse8387i" => "CCTGCA|GG"                   ,
    "sse8647i" => "AG|GWCCT"                    ,
    "sse9i" => "|AATT"                       ,
    "ssebi" => "AGG|CCT"                     ,
    "ssii" => "C|CGC"                       ,
    "sspd5i" => "GGTGA|"                      ,
    "sspdi" => "G|GCGCC"                     ,
    "sspi" => "AAT|ATT"                     ,
    "sste37i" => "CGAAGAC|"                    ,
    "ssti" => "GAGCT|C"                     ,
    "sth132i" => "CCCG|"                       ,
    "sth302ii" => "CC|GG"                       ,
    "stri" => "C|TCGAG"                     ,
    "stsi" => "GGATG|"                      ,
    "stui" => "AGG|CCT"                     ,
    "styd4i" => "|CCNGG"                      ,
    "styi" => "C|CWWGG"                     ,
    "swai" => "ATTT|AAAT"                   ,
    "taai" => "ACN|GT"                      ,
    "taii" => "ACGT|"                       ,
    "taqi" => "T|CGA"                       ,
    "taqii" => "GACCGA|"                     ,
    "tasi" => "|AATT"                       ,
    "tati" => "W|GTACW"                     ,
    "taui" => "GCSG|C"                      ,
    "tfii" => "G|AWTC"                      ,
    "tru1i" => "T|TAA"                       ,
    "tru9i" => "T|TAA"                       ,
    "tscai" => "CASTG|"                      ,
    "tsefi" => "|GTSAC"                      ,
    "tsei" => "G|CWGC"                      ,
    "tsoi" => "TARCCA|"                     ,
    "tsp45i" => "|GTSAC"                      ,
    "tsp4ci" => "ACN|GT"                      ,
    "tspdti" => "ATGAA|"                      ,
    "tspei" => "|AATT"                       ,
    "tspgwi" => "ACGGA|"                      ,
    "tspmi" => "C|CCGGG"                     ,
    "tspri" => "CASTG|"                      ,
    "tssi" => "GAGNNNCTC|GAGNNNCTC"         ,
    "tsti" => "CACN|NNNNNTCC"               ,
    "tsui" => "GCGAC|GCGAC"                 ,
    "tth111i" => "GACN|NNGTC"                  ,
    "tth111ii" => "CAARCA|"                     ,
    "ubaf11i" => "TCGTA|TCGTA"                 ,
    "ubaf12i" => "CTACNNNGTC|CTACNNNGTC"       ,
    "ubaf13i" => "GAGNNNNNNCTGG|GAGNNNNNNCTGG" ,
    "ubaf14i" => "CCANNNNNTCG|CCANNNNNTCG"     ,
    "ubaf9i" => "TACNNNNNRTGT|TACNNNNNRTGT"   ,
    "ubapi" => "CGAACG|CGAACG"               ,
    "ucomsi" => "|GAGCTC"                     ,
    "unbi" => "|GGNCC"                      ,
    "van91i" => "CCANNNN|NTGG"                ,
    "vha464i" => "C|TTAAG"                     ,
    "vnei" => "G|TGCAC"                     ,
    "vpak11ai" => "|GGWCC"                      ,
    "vpak11bi" => "G|GWCC"                      ,
    "vspi" => "AT|TAAT"                     ,
    "wvii" => "CACRAG|"                     ,
    "xagi" => "CCTNN|NNNAGG"                ,
    "xapi" => "R|AATTY"                     ,
    "xbai" => "T|CTAGA"                     ,
    "xcei" => "RCATG|Y"                     ,
    "xcmi" => "CCANNNNN|NNNNTGG"            ,
    "xhoi" => "C|TCGAG"                     ,
    "xhoii" => "R|GATCY"                     ,
    "xmai" => "C|CCGGG"                     ,
    "xmaiii" => "C|GGCCG"                     ,
    "xmaji" => "C|CTAGG"                     ,
    "xmii" => "GT|MKAC"                     ,
    "xmni" => "GAANN|NNTTC"                 ,
    "xspi" => "C|TAG"                       ,
    "ykri" => "C|"                          ,
    "zrai" => "GAC|GTC"                     ,
    "zrmi" => "AGT|ACT"                     ,
    "zsp2i" => "ATGCA|T"                     )
