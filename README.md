# MOSCA_webapp

1. *After Flask and MySQL installation here is what you should do:*
  - Open virtual environment with: 
    - $ . venv/bin/activate
  - Open MySQL shell and run following command: $ source <'absolute path for 'mosca_tables.sql''>. This will setup de database for the app, and create a test user for u, whose credetials are: Username: teste, Password:teste ;
  - Go back to linux, then change your directory for the folder that hold de file, *mosca_app.py*;
  - Run the following commands:
    - $ export FLASK_APP=mosca_app.py
    - $ flask run
    - Open the link generated by the last command, and the application should be good to go !
    
2. *If you intend to share the application with other users, the library local tunnel can be used to create a sharable URL that will work as long as the flask instance is running.
  -Install with :
    - $ npm install -g localtunnel
  - Open another terminal, and run local tunnel with the following command:
    - $ lt --port 5000
  
