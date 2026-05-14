#!/usr/bin/env python3
import sys
import yaml

file_path = sys.argv[1]
key = sys.argv[2]

with open(file_path) as f:
    config = yaml.safe_load(f)

value = config.get(key)
if value is None:
    sys.exit(f"ERRO: chave '{key}' não encontrada em {file_path}")
print(value)

