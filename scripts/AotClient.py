#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:10:50 2018

@author: guru

Array of Things API Query
"""

from aot_client import AotClient

client = AotClient()
projects = client.list_projects()
for page in projects:
  for proj in page.data:
    print(f'{proj["name"]} is available at /api/projects/{proj["slug"]}')