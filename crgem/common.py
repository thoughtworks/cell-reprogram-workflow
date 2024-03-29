"common.py"

import requests
import sys
import logging
import os
import shutil

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# artefacts_path = os.path.abspath(sys.argv[4])

def http_request(url, **kwargs):
    response = requests.get(url, **kwargs)
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
    return response

def http_post(url, **kwargs):
    return requests.post(url, **kwargs)

def clean_up(path):
    try:
        shutil.rmtree(path)
        logger.debug("Deleted alreday present subdirectory")
    except:
        logger.debug("No directory found to delete.")





