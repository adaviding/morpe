%================================================================
%function NTRIPS = Mcl_ForceMonotonic(OUTPUT_MONO, MEM_DY, Y, [IMETHOD])
% This function performs monotonic regression on Y.
% Warning:  The values of OUTPUT_MONO and MEM_DY will be changed inside this function.  Any data in these vectors will be destroyed.
%	This function outputs OUTPUT_MONO (a monotonic increasing function approximating Y).
%----------------------------------------------------------------
% OUTPUT_MONO:  (double) Outputs a non-decreasing 1-d vector the same size as Y and is a monotonic function fit to Y.
%	Guaranteed:
%		mean(OUTPUT_MONO) == mean(Y)
%
% MEM_DY:  (double) Allocated memory.  An array with one less element than Y.
%	
% Y:	(double) A 1-d vector.  This is the signal that will be forced monotonic (on output).  The input Y is not changed.
%
% IMETHOD:  An integer specifying the method to be used.  This argument is optional.  If no argument is supplied, default is 1.
%	(0):  A non-decreasing function is returned.  The function may may be flat in regions.
%	(1):  The non-decreasing function is interpolated through it's flat spots.
%			Each flat spot receives a proportion of positive derivative energy from immediately adjacent spots which are not flat.
%			The amount of energy received by the flat spots is proportional to the length (number of samples) of each flat spot.
%			If the non-decreasing function is flat everywhere (i.e. totally non-increasing), then a flat function will be returned.
%----------------------------------------------------------------
% NTRIPS:  The number of trips through a refining loop.
%================================================================