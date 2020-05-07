from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, SelectField, widgets, SelectMultipleField
from wtforms.validators import DataRequired
import pickle

class MultiCheckboxField(SelectMultipleField):
    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()

class SelectionForm(FlaskForm):
    geoaccession = StringField('GEOAccession', validators=[DataRequired()])

    ## load the chip information
    with open('app/website.pickle', 'rb') as f:
        data=pickle.load(f)

    chipType = SelectField(u'ChipType', choices=data )

    viz_options = MultiCheckboxField(u'Visualizations', choices=[('het', 'Heterozygosity'), ('pca', 'PCA'), ('circos', 'Circos'), ('sfs', 'Site Frequency Spectrum'), ('top', 'Heat Topology')])

    output = StringField('Output Label')

    samples = BooleanField("Click here if you want use all the samples in the study (or later you can specify specific samples)")

    submit = SubmitField('Done')

##class pickSamples(FlaskForm):
##    samples = MultiCheckboxField('Samples', choices=files)
##
##    submit = SubmitField('Done')