context("basic summaries")

test_that("basic summaries give correct numbers for iron data", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    expect_equal(n_ind(iron), 284)
    expect_equal(n_ind_geno(iron), 284)
    expect_equal(n_ind_pheno(iron), 284)
    expect_equal(n_ind_gnp(iron), 284)
    expect_equal(n_chr(iron), 20)
    expect_equal(chr_names(iron), paste0(c(1:19, "X")))
    expect_equal(tot_mar(iron), 66)
    expect_equal(n_mar(iron), c("1"=3, "2"=5, "3"=2, "4"=2, "5"=2, "6"=2, "7"=7, "8"=8,
                                "9"=5, "10"=2, "11"=7, "12"=2, "13"=2, "14"=2, "15"=2, "16"=5,
                                "17"=2, "18"=2, "19"=2, "X"=2))
    expect_equal(marker_names(iron), c("D1Mit18", "D1Mit80", "D1Mit17", "D2Mit379", "D2Mit75", "D2Mit17",
                                       "D2Mit304", "D2Mit48", "D3Mit22", "D3Mit18", "D4Mit2", "D4Mit352",
                                       "D5Mit11", "D5Mit30", "D6Mit104", "D6Mit15", "D7Mit74", "D7Mit25",
                                       "D7Nds5", "D7mit30", "D7Mit31", "D7Mit17", "D7Mit71", "D8Mit124",
                                       "D8Mit4", "D8Mit195", "D8Mit31", "D8Mit294", "D8Mit40", "D8Mit120",
                                       "D8Mit36", "D9Mit42", "D9Mit31", "D9Mit10", "D9Mit182", "D9Mit17",
                                       "D10Mit61", "D10Mit70", "D11Mit20", "D11Mit4", "D11Mit36", "D11Mit41",
                                       "D11Mit288", "D8Mit18", "D11Mit101", "D12Mit88", "D12Mit134",
                                       "D13Mit10", "D13Mit51", "D14Mit54", "D14Mit195", "D15Mit22",
                                       "D15Mit159", "D16Mit131", "D16Mit4", "D16Mit30", "D16Mit19",
                                       "D16Mit70", "D17Mit46", "D17Mit93", "D18Mit20", "D18Mit186",
                                       "D19Mit68", "D19Mit37", "DXMit16", "DXMit186"))
    expect_equal(n_pheno(iron), 2)
    expect_equal(pheno_names(iron), c("liver", "spleen"))
    expect_equal(n_covar(iron), 2)
    expect_equal(covar_names(iron), c("sex", "cross_direction"))
    expect_equal(n_phenocovar(iron), 1)
    expect_equal(phenocovar_names(iron), "tissue")
    expect_equal(ind_ids(iron), paste(1:284))
    expect_equal(ind_ids_geno(iron), paste(1:284))
    expect_equal(ind_ids_pheno(iron), paste(1:284))
    expect_equal(ind_ids_gnp(iron), paste(1:284))

    # run summary
    z <- summary(iron)

})

test_that("basic summaries give correct numbers for grav2 data", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    expect_equal(n_ind(grav2), 162)
    expect_equal(n_ind_geno(grav2), 162)
    expect_equal(n_ind_pheno(grav2), 162)
    expect_equal(n_ind_gnp(grav2), 162)
    expect_equal(n_chr(grav2), 5)
    expect_equal(chr_names(grav2), paste(1:5))
    expect_equal(tot_mar(grav2), 234)
    expect_equal(n_mar(grav2), c("1"=26, "2"=42, "3"=64, "4"=35, "5"=67))
    expect_equal(n_pheno(grav2), 241)
    expect_equal(pheno_names(grav2), paste0("T", seq(0, 480, by=2)))
    expect_equal(n_covar(grav2), 0)
    expect_equal(covar_names(grav2), NULL)
    expect_equal(n_phenocovar(grav2), 1)
    expect_equal(phenocovar_names(grav2), "time (hrs)")
    expect_equal(ind_ids(grav2), paste(1:162))
    expect_equal(ind_ids_geno(grav2), paste(1:162))
    expect_equal(ind_ids_pheno(grav2), paste(1:162))
    expect_equal(ind_ids_gnp(grav2), paste(1:162))

    expect_equal(marker_names(grav2), c("PVV4", "AXR-1", "HH.335C-Col/PhyA", "EC.480C", "EC.66C", "GD.86L",
                                        "CH.160L-Col", "CC.98L-Col/101C", "AD.121C", "AD.106L-Col", "GB.112L",
                                        "GD.97L", "EG.113L/115C", "CD.89C", "BF.206L-Col", "CH.200C",
                                        "EC.88C", "GD.160C", "HH.375L", "CH.215L", "BF.116C", "GH.157L-Col",
                                        "CC.318C", "CD.173L/175C-Col", "GH.127L-Col/ADH", "HH.360L-Col",
                                        "AD.156C", "BF.325L", "GH.580L", "DF.225L", "AD.77L", "CH.266C",
                                        "CH.610C", "HH.258C", "BH.145C", "BF.226C/BH.58L", "FD.226C",
                                        "GD.145C", "GH.94L", "BF.82C", "GD.465C", "FD.306L", "EC.495C-Col",
                                        "BH.460L", "FD.81L", "BF.105C", "CH.284C", "FD.222L-Col", "CD.245L",
                                        "EG.66L", "CH.65C", "CH.1500C", "BF.221L", "FD.85C", "GB.150L-Col",
                                        "FD.150C", "GD.460L-Col", "CC.332C", "Erecta", "CH.145L-Col/150C",
                                        "AD.191L-Col", "BH.195L-Col", "GD.298C", "GH.247L", "BH.120L-Col",
                                        "DF.140C", "EG.357C/359L-Col", "EC.235L-Col/247C", "DF.77C",
                                        "GB.120C-Col", "GD.248C-Col/249L", "EG.75L", "CH.322C", "FD.111L-Col/136C",
                                        "CC.266L", "BF.270L-Col/271C", "CC.110L/127C", "BH.88C", "EC.58C",
                                        "GH.390L", "GH.321L/323C-Col", "GH.226C/227L-Col", "HH.158L",
                                        "EC.83C/84L", "GD.318C/320L", "AD.92L", "HH.410C", "BF.148C",
                                        "GB.210L", "CD.800C", "HH.242C", "BF.134C-Col", "BH.225C-Col",
                                        "EG.122C", "BF.307L", "AD.427L", "BF.239L", "GH.411C/413L-Col",
                                        "GD.113C-Col", "BH.323C-Col", "GD.136C-Col", "GD.207C-Col", "GB.80C-Col",
                                        "AD.108L", "FD.97C-Col", "HH.440L", "GD.360L", "DF.303C", "DF.250L",
                                        "GD.174C-Col", "DF.76L", "BF.128C", "HH.117C", "HH.102C", "GD.296C-Col",
                                        "DF.328C", "GH.172C", "FD.98C", "CD.87L-Col", "CC.149L-Col",
                                        "AD.182C", "DF.65L-Col", "GD.106C", "GH.58C-Col", "HH.171C-Col/173L",
                                        "AD.179C", "AD.495L-Col", "AD.112L-Col", "BH.109L-Col", "GB.97L-Col/99C",
                                        "BH.285C", "HH.90L-Col", "ANL2", "GH.250C", "CH.169C", "FD.166L",
                                        "CD.320C", "BF.740L", "CC.288L", "GH.165L", "HH.161L", "CD.730C",
                                        "BF.151L", "AD.383C", "GH.266L", "CC.388L", "EC.306L", "GH.91C",
                                        "BH.92L-Col", "FD.154L", "GH.64C", "CH.238C", "CD.84C-Col/85L",
                                        "GH.263C-Col", "SC5", "g4539", "DF.108L-Col", "CD.329C-Col",
                                        "AD.307C", "AD.115L-Col", "FD.167L-Col", "CH.70L/71C-Col", "HH.159C-Col",
                                        "GH.433L-Col", "GB.490C", "GB.750C", "BH.342C/347L-Col", "FD.207L",
                                        "CH.690C", "AD.292L", "BH.144L", "FD.239L-Col", "EC.198L-Col",
                                        "BH.180C", "BH.325L", "BH.107L-Col", "BF.269C", "AD.114C-Col",
                                        "DF.231C", "CC.400L-Col", "DF.184L-Col", "EG.117L", "BF.164C",
                                        "GH.473C", "GH.117C", "GH.121L-Col", "DF.154C", "AD.129L-Col",
                                        "HH.480C", "BH.96L-Col", "CC.188L", "EC.395C", "CC.322C", "GH.190C",
                                        "GD.239L-Col", "CH.60C", "HH.225C-Col", "EG.410C", "DF.300C",
                                        "GB.235C-Col", "GB.59C", "CD.179L", "CH.88L", "EC.96L", "CD.116L",
                                        "CC.216C", "EC.151L", "CC.277L-Col", "DFR", "AD.254C", "AD.75C-Col",
                                        "GB.248C", "GD.350L-Col", "CD.160L", "BH.81L-Col", "BF.81C",
                                        "CC.540C-Col", "GB.223C", "DF.450C", "HH.445L-Col", "GD.118C",
                                        "CC.262C", "FD.345C", "GB.102L-Col/105C", "HH.143C", "BF.168L-Col",
                                        "DF.119L", "CH.331L-Col", "CH.124C", "GD.222C-Col", "g2368",
                                        "FD.188C", "EG.205L", "HH.122C/120L"))

    # run summary
    z <- summary(grav2)

})
