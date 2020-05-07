from flask import render_template, flash, redirect
from app import app
from app.forms import SelectionForm
import pipeline.gp_pipeline as pp
#import pipeline.visualization as pv

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='Home')

@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    form = SelectionForm()
    if form.validate_on_submit():
        flash("Data Submitted with GEO {} and ChipType {} for the graphs: {}".format(form.geoaccession.data, form.chipType.data, form.viz_options.data))
        flash('Use all the samples? {}'.format(form.samples.data))
        files = pp.main(form.geoaccession.data, form.chipType.data, form.samples.data, form.output.data)
        return redirect('/data')
    return render_template('analysis.html', title='Data Analysis', form=form)


@app.route('/data')
def data():
##    if form.samples.data == True:
##        return render_template('data.html', title='Data Running')
##    else:
##        return render_template('selectSamples.html', title='Select Samples')
    return render_template('data.html', title='Data Running')

@app.route('/visualization')
def visualization():
    return render_template('visualization.html', title='Visualization')

@app.route('/documentation')
def documentation():
    return render_template('documentation.html', title='Documentation')

@app.route('/faqs')
def faqs():
    return render_template('faqs.html', title='FAQs')

@app.route('/about')
def about():
    return render_template('about.html', title='About')