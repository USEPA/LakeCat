from datetime import datetime as dt
import json
import pandas as pd

title = 'The LakeCat Dataset: Accumulated Attributes for NHDPlusV2 (Version 2.1) Catchments for the Conterminous United States: '

KEYWORDS =['inlandwaters', 'ecosystem', 'environment', 'monitoring', 'natural resources',
   'surface water', 'modeling', 'united states of america', 'united states',
   'usa', 'alabama', 'arizona', 'arkansas', 'california', 'colorado', 'connecticut',
   'delaware', 'district of columbia', 'florida', 'georgia', 'idaho', 'illinois',
   'indiana', 'iowa', 'kansas', 'kentucky', 'louisiana', 'maine', 'maryland',
   'massachusetts', 'michigan', 'minnesota', 'mississippi', 'missouri', 'montana',
   'nebraska', 'nevada', 'new hampshire', 'new jersey', 'new mexico', 'new york',
   'north carolina', 'north dakota', 'ohio', 'oklahoma', 'oregon', 'pennsylvania',
   'rhode island', 'south carolina', 'south dakota', 'tennessee', 'texas', 'utah',
   'vermont', 'virginia', 'washington', 'west virginia', 'wisconsin', 'wyoming']

organization = 'U.S. Environmental Protection Agency, Office of Research and Development (ORD), Center for Public Health and Environmental Assessment (CPHEA), Pacific Ecological Systems Division (PESD), '
temporal = '2015/2030'

LakeCat = ('LakeCat currently contains over 300 metrics that include local catchment (Cat),'
             ' watershed (Ws), and special metrics. See Geospatial Framework and Terms in the'
             ' ReadMe for definitions of the terms \u2018catchment\u2019 and \u2018watershed\u2019'
             ' as used with the LakeCat Dataset. An additional metric, inStreamCat, indicates whether'
             ' the variable was pulled from the StreamCat Dataset or calculated with a geospatial'
             ' framework that was developed for LakeCat.\n\nThese metrics are available for 378,088'
             ' lakes and their associated catchments across the conterminous US. LakeCat metrics'
             ' represent both natural (e.g., soils and geology) and anthropogenic (e.g, urban areas'
             ' and agriculture) landscape information.')


describeBy = 'https://www.epa.gov/national-aquatic-resource-surveys/lakecat-metrics-and-definitions' #
access = 'https://www.epa.gov/national-aquatic-resource-surveys/lakecat-dataset'
gaft =  'https://gaftp.epa.gov/epadatacommons/ORD/NHDPlusLandscapeAttributes/LakeCat/FinalTables/' 


tbl = pd.read_csv(
    "O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/Metadata"
    "/submit_metadata_LakeCat.csv" #replace
)

records = []

for _, row in tbl.iterrows():
    # break
    records.append(
        {'@type': "dcat:Dataset",
         'title': title + row.title,                                           #made title variable
         'description': row.description,
         'keyword': KEYWORDS,
         'modified': dt.now().strftime('%Y_%m_%d'),
         'publisher':
             {'@type': 'org:Organization',
              'name': organization                                             # made organization variable
         },
         'contactPoint': {
         '@type': 'vcard:Contact',
         'fn': organization + row.contact,                                     # made 'contact' variable and fill "Marc Weber"
         'hasEmail': ('mailto:' + row.contact_Email)                           # made 'contact_Email' variable and fill weber.marc@epa.gov
         },
         'identifier': row.uuid,                                               # got rid of f"{{{row.uuid}}}" in favor of row.uuid
         'accessLevel': "public",
         'bureauCode': ['020:00'],
         'programCode': ['020:072'],
         'license': "https://edg.epa.gov/EPA_Data_License.html",
         'rights': "public (Data asset is or could be made publicly available to all without restrictions)",
         'spatial': '-125.0,24.5,-66.5,49.5',
         'temporal': temporal,
         'distribution': [
            {'@type': 'dcat:Distribution',
             'accessURL': access,                                              #
             'title': 'LakeCat Dataset',
             'description': LakeCat,                                           #made LakeCat a variable to make cleaner
			 'format': 'API',                                                  #changed to API
             'describedBy': describeBy,                                        # made URL variable
             'describedByType': 'text/html'
            },
            {'@type': 'dcat:Distribution',
             'downloadURL': gaft,                                              #
             'format': 'Comma-Separated Values (.csv)',
             'title': row.final_table_name,
             'description': row.description,
             "mediaType": "text/csv",
             "describedBy": describeBy,
          "describedByType": "text/html"
             }
          ],
         'accrualPeriodicity': "irregular",
         'dataQuality': True,
         'describedBy': describeBy,
         'describedByType': "text/html",
         'issued': "2015-04-23",
         'accrualPeriodicity': "R/P3Y",
         'language': ['en-US'],
         'landingPage': access,                                                #
         'references': [
             (gaft),                                                             #
             ('https://edg.epa.gov/metadata/catalog/search/resource/details.page?uuid=' + row.uuid)  ###removed f'{{{row.uuid}}}
         ],
         'theme': ['environment'],
        }
    )

blob = {
      "@context": 'https://project-open-data.cio.gov/v1.1/schema/catalog.jsonld', #
      "@id": "https://edg.epa.gov/data/Public/ORD/NHEERL/WED/NonGeo/LakeCat.json", #
      "@type": "dcat:Catalog",
      "conformsTo": "https://project-open-data.cio.gov/v1.1/schema", #check validation in open
      "describedBy": "https://project-open-data.cio.gov/v1.1/schema/catalog.json",
      "dataset": records,
}
with open(
        "O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/Metadata/copmlete_metadata"
        f"/LakeCat_metadata_{dt.now().strftime('%m_%d_%Y')}.json", #replace 
        "w"
    ) as fifi:
    json.dump(blob, fifi, sort_keys=False, indent=4,)          

                                                                   


