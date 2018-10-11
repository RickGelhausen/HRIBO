rule xtailreport:
    input:
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
    output:
        fcplot=report("figures/fc_{contrast}.jpg", caption="../report/xtail_fc.rst", category="Regulation - Xtail"),
        rplot=report("figures/r_{contrast}.jpg", caption="../report/xtail_r.rst", category="Regulation - Xtail")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 150 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; convert -density 150 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 150 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; ")

rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf" 
    output:
        report("figures/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Novel ORFs - Ribotish")
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
        metagene=report("figures/{condition}-{replicate}_metagene.jpg", caption="../report/reparation_metagene.rst", category="Novel ORFs - Reparation"),
        roc=report("figures/{condition}-{replicate}_roc.jpg", caption="../report/reparation_roc.rst", category="Novel ORFs - Reparation"),
        psite=report("figures/{condition}-{replicate}_psite.jpg", caption="../report/reparation_psite.rst", category="Novel ORFs - Reparation"),
        scurve=report("figures/{condition}-{replicate}_scurve.jpg", caption="../report/reparation_scurve.rst", category="Novel ORFs - Reparation")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input.metagene} -quality 100 -flatten -sharpen 0x1.0 {output.metagene}; convert -density 150 -trim {input.roc} -quality 100 -flatten -sharpen 0x1.0 {output.roc}; convert -density 150 -trim {input.psite} -quality 100 -flatten -sharpen 0x1.0 {output.psite}; convert -density 150 -trim {input.scurve} -quality 100 -flatten -sharpen 0x1.0 {output.scurve};")
