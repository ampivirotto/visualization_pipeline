import os
from igvjs import app
from flask_ngrok import run_with_ngrok

run_with_ngrok(app)
#port = int(os.environ.get('PORT', 5000))
#host = '127.0.0.1' if app.config['DEBUG'] else '0.0.0.0'
app.run()
