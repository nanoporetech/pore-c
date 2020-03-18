Data Model
===========


Fragment Alignments
--------------------

This table contains the Pore-C alignments, whether they've been filtered out and which fragments they've been assigned to.


.. jsonschema:: data_model.json#/definitions/PoreCRecord


Contacts and Concatemers
--------------------------

These tables are derived from the Fragment Alignment records. The PoreCContactRecord represents a single pairwise contact and the PoreCConcatemerRecord contains aggregate statistics for all of the contacts for a read.


.. jsonschema:: data_model.json#/definitions/PoreCContactRecord
.. jsonschema:: data_model.json#/definitions/PoreCConcatemerRecord
