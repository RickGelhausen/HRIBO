rule contrastInput:
    output:
        "contrasts/{contrast}"
    run:
        if not os.path.exists("contrasts"):
            os.makedirs("contrasts")
        for f in getContrast(wildcards):
            print(f)
            open((f), 'a').close()
