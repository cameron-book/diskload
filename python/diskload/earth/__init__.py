from collections import namedtuple

Model = namedtuple("Model", "radius gravity")
default = Model( radius=6371.0, gravity=9.8046961 )
