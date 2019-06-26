import argparse
from passlib.hash import sha256_crypt
from flask_mysqldb import MySQL
import mysql.connector



parser = argparse.ArgumentParser(description="""Script for adding new users to MOSCA's webapp""")
parser.add_argument("-n","--name", type=str,
                    help="Name of new user")
parser.add_argument("-e","--email", type=str,
                    help="email of new user")
parser.add_argument("-u","--username", type=str,
                    help="Username of new user")
parser.add_argument("-p","--password", type=str,
                    help="Password for new user")

parser.add_argument("-c","--command" type=str, help="action to perform", default = "add")

args = parser.parse_args()

cnx = mysql.connector.connect(
  host="localhost",
  user="root",
  passwd="",
  database="myproject"
)

def insert_new_user(name, username, password, email=''):


    #cnx = mysql.connector.connect(user='root',database='myproject')
    cur = cnx.cursor(dictionary=True)

    password = sha256_crypt.encrypt(str(password))

    cur.execute('INSERT INTO users(name, email, username, password) VALUES(%s, %s, %s, %s)',(name, email, username, password))

    # COMIT to db
    cnx.commit()
    #close connection
    cur.close()

def remove_user(username):
  cur = cnx.cursor(dictionary=True)
  cur.execute('DELETE FROM users WHERE username=%s',[username])
  # COMIT to db
  cnx.commit()
  #close connection
  cur.close()

  
if args.command == 'add':
  insert_new_user(args.name, args.username, args.password,args.email)
elif args.command=='remove':
  remove_user(args.username)
else:
  print('Command is either add or remove')
