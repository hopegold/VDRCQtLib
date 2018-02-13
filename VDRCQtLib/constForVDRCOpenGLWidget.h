#pragma once

#include "Color3f.h"

const int SPHERE_RESOLUTION = 100;

const float VEDGE_THICKNESS = 1.0f;

enum  PICK_MODE
{
	NON_PICKING,
	PICK_VVERTEX
};

const int NUM_PICKING_CLASS = 10;
const int CLASS_VVERTEX = 1;

const float SELECTION_BOX_SIZE = 10.0;
const int SELECTION_BUFFER_SIZE = 100;

const float VVERTEX_BALL_RADIUS = 1;

const float STIPPLE_FACTOR = 4;
const GLushort STIPPLE_PATTERN = 0xAAAA;

const float CONE_HEIGHT = 1;
const float CONE_BASE_RADIUS = 0.5;
const float SPEED_RATIO_FOR_DRAWING_CONE = 100;