"common.py"

import requests
import sys
import datetime

def http_request(url, **kwargs):
    response = requests.get(url, **kwargs);
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
    return response

def http_post(url, **kwargs):
    return requests.post(url, **kwargs)


# def unique_directory_name():
#     return "Results-"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
