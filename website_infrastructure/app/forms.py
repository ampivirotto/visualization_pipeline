from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired

class SelectionForm(FlaskForm):
    geoaccession = StringField('GEOAccession', validators=[DataRequired()])

    chipType = SelectField('ChipType', choices =[("human_omni_express", "HumanOmniExpress-24-v1-0-B"), ("dogHD","CanineHD"), ("inf_omni_zhonghua", "InfiniumOmniZhongHua"),
                                                    ("human_omni", "HumanOmni5-4")])

    #location = StringField('Local Location', validator=[DataRequired()])

    submit = SubmitField('Done')