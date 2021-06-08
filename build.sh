#!/bin/bash

echo "--------------------Uninistall any old versions...-----------------"
pip3 uninstall -y pycheron

echo "--------------------Building Pycheron----------------------"
poetry build

echo "--------------------Installing Pycheron--------------------"
poetry run pip3 install --upgrade pip
poetry run pip3 install numpy==1.19.1
poetry run pip3 install llvmlite==0.34.0
poetry run pip3 install ./dist/pycheron-3.0.0-py3-none-any.whl
