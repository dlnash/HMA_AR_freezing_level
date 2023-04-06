#!/bin/bash
######################################################################
# Filename:    dryad_submission.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to get token API access and upload data to dryad repo
#
######################################################################

# Input parameters
domain_name="datadryad.org" # domain name
# application_id="dhFE0xvzWOA6zgSpGTuapASUx0h4hKFyk5eNvRgqb_M"
# secret="YBZm-aeQJRHwpL3u68lQkv7RpXHql7kaMRrL5tF38Tc"
access_token="CRVePQYpSiLPr0FK0gIkBBTY0_SC6NQIw3DFCfkpOJQ"
file_path="/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/"
filename="README.md"
doi="10.25349/D9D61Q"
content_type="text/plain" #application/zip

encoded_doi="doi%10.25349%D9D61Q"

encoded_filename='README.md' | base64

# python -c 'import sys,pathlib; print(pathlib.Path(sys.argv[1]).resolve().as_uri())' "https://doi.org/10.25349/D9D61Q"

# # get token with curl
# curl -X POST https://${domain_name}/oauth/token -d "client_id=${application_id}&client_secret=${secret}&grant_type=client_credentials" -H "Content-Type: application/x-www-form-urlencoded;charset=UTF-8"

# ## returns
# {"access_token":"CRVePQYpSiLPr0FK0gIkBBTY0_SC6NQIw3DFCfkpOJQ","token_type":"Bearer","expires_in":36000,"scope":"all","created_at":1680712303}


## test token
# curl -i -H "Accept: application/json" -H "Content-Type: application/json" -H "Authorization: Bearer ${access_token}" -X GET https://${domain_name}/api/v2/test

## test upload README.md
# curl --data-binary "@</path/to/my/file>" -i -X PUT "https://<domain-name>/api/v2/datasets/<encoded-doi>/files/<encoded-file-name>" -H "Authorization: Bearer <token>" -H "Content-Type: <mime-type>" -H "Accept: application/json"


curl --data-binary "@${file_path}${filename}" -i -X PUT "https://${domain_name}/api/v2/datasets/${encoded_doi}/files/${encoded_filename}" -H "Authorization: Bearer ${access_token}" -H "Content-Type: ${content_type}" -H "Accept: application/json"

