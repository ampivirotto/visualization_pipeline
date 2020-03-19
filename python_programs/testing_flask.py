#######################  IMPORT #############################################
from flask_wtf import FlaskForm
from flask import Flask
from wtforms import StringField, BooleanField, SubmitField
from wtforms.validators import DataRequired
app = Flask(__name__)


##############################  FUNCTIONS #####################################


class SelectionForm(FlaskForm):
    geoaccession = StringField('GEOAccession', validator=[DataRequired()])
    location = StringField('Local Location', validator=[DataRequired()])
    submit = SubmitField('Done')



################################  MAIN #####################################
if __name__ == '__main__':
    app.run()