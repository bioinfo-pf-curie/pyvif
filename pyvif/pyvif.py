# -*- coding: utf-8 -*-

"""Main module."""

from jinja2 import Environment, PackageLoader

from pyvif.utils import include_file, template_data


def main():
    env = Environment(loader=PackageLoader('pyvif', 'templates'))
    env.globals['include_file'] = include_file
    env.globals['template_data'] = template_data
    template = env.get_template('main.html')
    report_output = template.render()

    with open('test.html', 'w') as fout:
        print(report_output, file=fout)
