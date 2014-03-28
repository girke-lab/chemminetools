#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import json
import sys
import os

parser = \
    argparse.ArgumentParser(description='Generate cytoscape plot'
                            )
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
parser.add_argument('-l', '--layout', help='layout type',required=True)
args = vars(parser.parse_args())

def main():
    # parse input
    fa = sys.stdin.read()
    fa = fa.rstrip()

    outputData = {
            u'maxZoom': 2, 
            u'minZoom': 0.5, 
            u'style': {
                u'1': 
                    {u'properties': 
                        [{u'name': u'text-outline-color', u'value': u'#000'}, {u'name': u'background-color', u'value': u'#000'}, {u'name': u'target-arrow-color', u'value': u'#000'}, {u'name': u'line-color', u'value': u'#000'}], 
                        u'selector': u':selected'}, 
                    u'0': {u'properties': [{u'name': u'text-valign', u'value': u'center'}, {u'name': u'color', u'value': u'#fff'}, {u'name': u'content', u'value': u'data(name)'}, {u'name': u'text-outline-color', u'value': u'#888'}, {u'name': u'text-outline-width', u'value': 3}, {u'name': u'font-family', u'value': u'helvetica'}, {u'name': u'font-size', u'value': 14}, {u'name': u'border-color', u'value': u'#fff'}, {u'name': u'height', u'value': u'mapData(height, 0, 200, 10, 45)'}, {u'name': u'width', u'value': u'mapData(weight, 30, 80, 20, 50)'}], u'selector': u'node'}, 
                    u'length': 3, u'2': {u'properties': [{u'name': u'width', u'value': 2}, {u'name': u'target-arrow-shape', u'value': u'triangle'}], u'selector': u'edge'}
            }, 
            u'showOverlay': False
            }

    outputData['elements'] = json.loads(fa)
    outputData['layout'] = {
        u'name': u'arbor'
    }
    outputData['style'] = [ 
        {'selector': 'node',
            'css': {
            'content': 'data(name)'
            }}
    ]

    # write options to output file
    of = open(args['outfile'], 'w')
    of.write(json.dumps(outputData))
    of.close()

if __name__ == '__main__':
    main()
