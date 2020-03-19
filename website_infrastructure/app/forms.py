from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField
from wtforms.validators import DataRequired

class SelectionForm(FlaskForm):
    geoaccession = StringField('GEOAccession', validator=[DataRequired()])
    #location = StringField('Local Location', validator=[DataRequired()])
    submit = SubmitField('Done')