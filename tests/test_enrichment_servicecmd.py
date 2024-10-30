#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `enrichment_servicecmd` package."""

import os
import tempfile
import shutil

import unittest
from enrichment_service import enrichment_servicecmd


class TestEnrichmentServiceCommand(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_genes_from_data_with_list(self):
        self.assertEqual([], enrichment_servicecmd.get_genes_from_data([]))
        self.assertEqual(['a', 'b'],
                         enrichment_servicecmd.get_genes_from_data([' a ', 'b']))

    def test_get_genes_from_data_with_str(self):
        self.assertEqual([], enrichment_servicecmd.get_genes_from_data(''))
        self.assertEqual(['a', 'b', 'c', 'dss'],
                         enrichment_servicecmd.get_genes_from_data(' a b,c , dss  '))
