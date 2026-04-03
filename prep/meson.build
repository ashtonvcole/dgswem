prep_src = files([
'adcprep.F',
'decomp.F',
'interp.F',
'machdep.F',
'metis.F',
'metis2.F',
#'picomsg.F', This file is still here and included in the makefile, but isn't used by the makefile? bizzare.
'pre_global.F',
'prep.F',
'presizes.F',
'read_global.F'
])

post_src = files([
  'adcpost.F',
  'compare.F',
  'diffmerge.F',
  'post.F',
  'post_global.F'])


metis_dep = dependency('metis', version: '>=4.0.1',
                       required: true,
                       fallback: ['metis', 'metis_dep'])
 

