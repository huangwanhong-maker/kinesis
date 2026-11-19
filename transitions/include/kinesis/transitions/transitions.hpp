// kinesis::transitions — critical transitions read from records of structured states.
//
// space.hpp       records, state spaces and their charts (identity, composition, histogram, covariance)
// linear.hpp      covariance, the lag-1 propagator, leading modes, Lyapunov equations
// rate.hpp        the restoring rate by four estimators
// nulls.hpp       the fitted null, its slowed alternative, and surrogates
// trend.hpp       trends in indicators with their operating characteristics, robustness and grid
// kinds.hpp       where the dominant eigenvalue goes
// states.hpp      wells, switching and Kramers's time
// likelihood.hpp  a drifting rate against a constant one, with a tipping time only when it is bounded
// indicators.hpp  the indicators with what they see and are blind to
// detect.hpp      an assessment with verdicts and reasons
// models.hpp      records with known answers
#pragma once

#include <kinesis/transitions/detect.hpp>
#include <kinesis/transitions/models.hpp>
