

def Uniprot2String(identifiers, species=9606):

    ##########################################################
    ## For a given list of proteins the script resolves them
    ## (if possible) to the best matching STRING identifier
    ## and prints out the mapping on screen in the TSV format
    ##
    ## Requires requests module:
    ## type "python -m pip install requests" in command line
    ## (win) or terminal (mac/linux) to install the module
    ###########################################################

    import requests ## python -m pip install requests

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    ##
    ## Set parameters
    ##

    params = {

        #"identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "identifiers" : "\r".join(identifiers), # your protein list
        #"species" : 9606, # species NCBI identifier 
        "species" : species, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "www.awesome_app.org" # your app name

    }

    ##
    ## Construct URL
    ##


    request_url = "/".join([string_api_url, output_format, method])

    ##
    ## Call STRING
    ##

    results = requests.post(request_url, data=params)

    ##
    ## Read and parse the results
    ##

    res = []
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, string_identifier = l[0], l[2]
        res.append((input_identifier, string_identifier))
        #print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
    return res


def FunctionalEnrichment(my_genes, background, species=9606):

    ##############################################################
    ## The following script retrieves and prints out
    ## significantly enriched (FDR < 1%) GO Processes
    ## for the given set of proteins. 
    ##
    ## Requires requests module:
    ## type "python -m pip install requests" in command line (win)
    ## or terminal (mac/linux) to install the module
    ##############################################################

    import requests ## python -m pip install requests 
    import json

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv"
    method = "enrichment"


    ##
    ## Construct the request
    ##

    request_url = "/".join([string_api_url, output_format, method])

    ##
    ## Set parameters
    ##

    # my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
    #              '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

    params = {

        "identifiers" : "%0d".join(my_genes), # your protein
        'background_string_identifiers': "%0d".join(background),
        # "species" : 7227, # species NCBI identifier 
        #"species" : species, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org" # your app name

    }

    ##
    ## Call STRING
    ##

    response = requests.post(request_url, data=params)

    ##
    ## Read and parse the results
    ##

    return response.text
    data = json.loads(response.text)
    for row in data:

        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = float(row["fdr"])
        description = row["description"]
        category = row["category"]

        #if category == "Process" and fdr < 0.01:
        if fdr < 0.01:

            ## print significant GO Process annotations

            print("\t".join([term, preferred_names, str(fdr), description, category]))