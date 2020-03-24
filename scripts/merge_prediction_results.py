# Used to create a common annotation to get common read count stats

def create_merged_annotation(args):
    if args.deepribo_path == "":
        
    else:

        

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create an annotation file from deepribo and reparation.')
    parser.add_argument("-d", "--deepribo_predictions", action="store", dest="deepribo_path", default="", help= "deepribo predictions gff file.")
    parser.add_argument("-r", "--reparation_predictions", action="store", dest="reparation_path", required=True, help= "reparation predictions gff file.")
    parser.add_argument("-o", "--prediction_annotation", action="store", dest="annotation_path", required=True, help= "annotation containing gff file.")
    args = parser.parse_args()

if __name__ == '__main__':
    main()
