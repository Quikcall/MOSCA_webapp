#!/bin/bash
mkdir MOSCA_app
cd MOSCA_app
python3 -m venv venv
. venv/bin/activate
pip install Flask
apt-get install default-libmysqlclient-dev
pip install flask_mysqldb
pip install wtforms
pip install passlib
pip install mysql-connector
mysql -u root -p -e 'source mosca_tables.sql'

