from bs4 import BeautifulSoup
import json

def main(homer_html, save_file):
    html_doc = open(homer_html, "r")
    soup = BeautifulSoup(html_doc, 'html.parser')

    motif_dict = dict()

    i=0

    for tag_row in soup.find_all("tr"):
        if i == 0:
            i+=1
            continue
        motif_info = tag_row.contents[1].contents
        motif_name = motif_info[0].rstrip("(1.000)")
        if len(motif_info) == 5:
            motif_dict[motif_name] = motif_info[-1].contents[0].rstrip(" info")
        else:
            motif_dict[motif_name] = ""
        i+=1  
    
    with open(save_file, "w") as outfile:
        json.dump(motif_dict, outfile, indent=4)

    return


if __name__ == "__main__":
    homer_html = "/data5/deepro/starrseq/main_library/10_enhancer_characterization/data/homer/homerResults.html"
    save_file = "/data5/deepro/starrseq/main_library/10_enhancer_characterization/data/homer/homerResults.json"
    main(homer_html, save_file)