/*
  genetic.c
  $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
  $Revision: 134 $
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "models.h"
#include "missing.h"
#include "recommend.h"
#include "genetic.h"

#define max(A, B) ((A > B)? A : B)

