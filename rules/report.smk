rule xtailreport:
    input:
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
    output:
        fcplot=report("figures/fc_{contrast}.jpg", caption="../report/xtail_fc.rst", category="Regulation"),
        rplot=report("figures/r_{contrast}.jpg", caption="../report/xtail_r.rst", category="Regulation")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 100 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 100 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; convert -density 100 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 100 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; ")

rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf"
    output:
        "figures/{condition}-{replicate}-qual.jpg"
        #report("figures/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Novel ORFs - Ribotish")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input} -quality 100 -flatten -sharpen 0x1.0 {output};")

rule reparationreport:
    input:
        metagene="reparation/{condition}-{replicate}/metagene_profile.pdf",
        roc="reparation/{condition}-{replicate}/PR_and_ROC_curve.pdf",
        psite="reparation/{condition}-{replicate}/p_site_offset.png",
        scurve="reparation/{condition}-{replicate}/S_Curve.pdf"
    output:
        metagene="figures/{condition}-{replicate}_metagene.jpg",
        roc="figures/{condition}-{replicate}_roc.jpg",
        psite="figures/{condition}-{replicate}_psite.png",
        scurve="figures/{condition}-{replicate}_scurve.jpg"
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input.metagene} -quality 100 -flatten -sharpen 0x1.0 {output.metagene}; convert -density 150 -trim {input.roc} -quality 100 -flatten -sharpen 0x1.0 {output.roc}; cp {input.psite} {output.psite}; convert -density 150 -trim {input.scurve} -quality 100 -flatten -sharpen 0x1.0 {output.scurve};")
