// The grammar of the language of dynamical systems, in BNF over bytes, read by morphe's reader
// (morphe/formal/bnf.hpp): the first rule is the start, names beginning with '_' are hidden, and whitespace is
// written out, since the grammar is over bytes.
//
// A script is a list of systems. A system is a declaration in parts, or a derivation from another system:
//
//   state, param, noise, unknown, input   the carrier S, its parameters, its sources of noise, what is not yet
//                                         known, and the signals it takes from outside
//   change                                the generators and how each acts: equations (a flow x' = f, a map, a
//                                         stochastic differential dX = a*dt + b*dW), an action, a Hamiltonian
//                                         with its jumps, rules, and events with a guard and resets
//   structure                             how changes compose into histories and what a history weighs: the
//                                         time, derivations, the reading of an action, the semiring of weights,
//                                         the stochastic calculus
//   constraint                            what every change respects: a background, an invariance, a conserved
//                                         quantity, an identity
//   difference                            what is told apart: the observations, or what states are the same up to
//   part, connect                         a system made of others: its parts, and what each input is given
//
// Choices (the idea of 2026-11-12, §9.1). Blocks are in braces. A product is written with '*': juxtaposition would
// be ambiguous with application. A line may continue before a binary operator, and no item begins with one, the
// left side of an equation being a name, a primed name or an application. The draft is ASCII only.
//
// A change to this text changes the reader (reader.hpp), the declaration (declaration.hpp) and the checks together.
#pragma once

#include <string_view>

namespace kinesis::language {

inline constexpr std::string_view grammar_text = R"BNF(
script    ::= _ ( system ( _sep system )* _ )? ;

_    ::= ( [ \t\n\r] | _c )* ;
_c   ::= "#" [^\n]* ;
_s   ::= [ \t]* ;
_sp  ::= [ \t]+ ;
_sep ::= _s _c? ( [\n\r;] _s _c? )+ ;

name   ::= [A-Za-z_] [A-Za-z0-9_]* ;
number ::= [0-9]+ ( "." [0-9]+ )? ( [eE] [+\-]? [0-9]+ )? ;
names  ::= name ( _s "," _s name )* ;

system     ::= "system" _sp name _s ( body | "=" _s derivation ) ;
body       ::= "{" _ ( _item ( _sep _item )* _ )? "}" ;
_item      ::= state | param | noise | unknown | input | part | change | structure | constraint | difference | connect ;
derivation ::= dop _sp name ;
dop        ::= "expand" | "reverse" | "lift" | "vary" | "legendre" | "unravel" ;

state     ::= "state" _sp decls _s ":" _s space ( _s "," _s qualifier )* ;
decls     ::= decl ( _s "," _s decl )* ;
decl      ::= name ( "(" _s names _s ")" )? ;
space     ::= wordspace | name ( "(" _s args _s ")" )? ( "^" exponent )? ;
wordspace ::= "words" _sp "over" _s set ;
exponent  ::= number | name ;
set       ::= "{" _s names _s "}" ;
qualifier ::= ranging | start ;
ranging   ::= name ( _sp name )* _sp "in" _sp domain ;
start     ::= "start" _sp word ;
domain    ::= interval | group | eapp | name ;
group     ::= "group" _sp name ;
interval  ::= "[" _s expr _s "," _s expr _s "]" ;

param    ::= "param" _sp binding ( _s "," _s binding )* ;
binding  ::= decl _s "=" _s expr ;
noise    ::= "noise" _sp target _s nkind _s expr ;
nkind    ::= ":" | "~" ;
target   ::= name ( "(" _s args _s ")" )? ;
unknown  ::= "unknown" _sp names _s ":" _s ukind ;
ukind    ::= "number" | "function" | "sparse" _sp expr ;
input    ::= "input" _sp ibinding ( _s "," _s ibinding )* ;
ibinding ::= name ( _s "=" _s expr )? ;

part    ::= "part" _sp names _s ":" _s name ( _sp "with" _sp binding ( _s "," _s binding )* )? ;
connect ::= "connect" _s "{" _ ( wire ( _sep wire )* _ )? "}" ;
wire    ::= qname _s "=" _s expr ;
qname   ::= name "." name ;

change      ::= "change" _s "{" _ ( _citem ( _sep _citem )* _ )? "}" ;
_citem      ::= action | hamiltonian | jump | event | rule | equation ;
action      ::= "action" _sp name _s "=" _s expr ;
hamiltonian ::= "hamiltonian" _sp name _s "=" _s expr ;
jump        ::= "jump" _sp name _s "=" _s expr ;
event       ::= "when" _sp expr _s relop _s expr _s ":" _s assign ( _s "," _s assign )* ;
relop       ::= ">=" | "<=" | ">" | "<" ;
assign      ::= name primes? _s "->" _s expr ;
rule        ::= name _s ":" _s word _s "->" _s word ( _s "@" _s expr )? ;
word        ::= name ( _sp name )* ;
equation    ::= lhs _s "=" _s expr ;
lhs         ::= name ( primes | "(" _s args _s ")" )? ;
primes      ::= "'"+ ;

structure   ::= "structure" _s "{" _ ( _sitem ( _sep _sitem )* _ )? "}" ;
_sitem      ::= time | derivations | histories | weights | calculus ;
time        ::= "time" _sp name _sp "in" _sp domain ;
derivations ::= "derivations" ;
histories   ::= "histories" _sp ( stationary | weight ) ;
stationary  ::= "stationary" _sp name ;
weight      ::= "weight" _sp expr ;
weights     ::= "weights" _sp wkind ;
wkind       ::= "boolean" | "probability" | "real" | "complex" | "operators" ;
calculus    ::= "calculus" _sp ckind ;
ckind       ::= "ito" | "stratonovich" ;

constraint ::= "constraint" _s "{" _ ( _kitem ( _sep _kitem )* _ )? "}" ;
_kitem     ::= background | invariance | conserved | equation ;
background ::= name ( _sp name )* _sp "on" _sp domain ;
invariance ::= "invariant" _sp "under" _sp expr ;
conserved  ::= "conserves" _sp expr ;

difference ::= "difference" _s "{" _ ( _ditem ( _sep _ditem )* _ )? "}" ;
_ditem     ::= observe | upto ;
observe    ::= "observe" _sp obs ( _s "," _s obs )* ;
obs        ::= ( name _s "=" _s )? expr ;
upto       ::= "up" _sp "to" _sp expr ;

expr    ::= _eadd ;
_eadd   ::= _emul ( _ addop _s _emul )* ;
addop   ::= "+" | "-" ;
_emul   ::= _eunary ( _ mulop _s _eunary )* ;
mulop   ::= "*" | "/" ;
_eunary ::= neg | _epow ;
neg     ::= "-" _s _eunary ;
_epow   ::= power | _eatom ;
power   ::= _eatom "^" _eunary ;
_eatom  ::= number | hole | eapp | qvar | evar | paren ;
qvar    ::= qname primes? ;
paren   ::= "(" _ expr _ ")" ;
eapp    ::= name "(" _s ( args _s )? ")" ;
evar    ::= name primes? ;
hole    ::= "?" name ( "(" _s args _s ")" )? ;
args    ::= expr ( _s "," _s expr )* ;
)BNF";

}  // namespace kinesis::language
