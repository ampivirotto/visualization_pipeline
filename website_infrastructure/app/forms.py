from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired

class SelectionForm(FlaskForm):
    geoaccession = StringField('GEOAccession', validators=[DataRequired()])

    chipType = SelectField('ChipType', choices =[("humanOmni", "HumanOmniExpress-24-v1-0-B"), ("dogHD","CanineHD")])

    #location = StringField('Local Location', validator=[DataRequired()])

    submit = SubmitField('Done')