from flask import render_template, flash, redirect
from app import app
from app.forms import SelectionForm

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='Home')

@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    form = SelectionForm()
    if form.validate_on_submit():
        flash('Data Submitted with GEO {} and ChipType {}'.format(form.geoaccession.data, form.chipType.data))
        return redirect(url_for('/data'))
    return render_template('analysis.html', title='Data Analysis', form=form)

@app.route('/data')
def data():
    return render_template('data.html', title='Data in Database')

@app.route('/visualization')
def visualization():
    return render_template('visualization.html', title='Visualization')

@app.route('/documentation')
def documentation():
    return render_template('documentation.html', title='Documentation')