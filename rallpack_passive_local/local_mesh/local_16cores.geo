//+
SetFactory("OpenCASCADE");
Cylinder(1) = {1.6, 1.6, 0, 10, 0, 0, 0.5, 2*Pi};
//+
Cylinder(2) = {1.7, -1.4, 0, 10, 0, 0, 0.5, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Volume{2}; 
}
//+
Translate {0, 0, 1} {
  Volume{2}; 
}
//+
Translate {0, 0, 1} {
  Volume{2}; 
}
//+
Translate {0, 1, 0} {
  Volume{2}; 
}
//+
Translate {0, 1, 0} {
  Volume{2}; 
}
//+
Translate {0, 1, 0} {
  Volume{2}; 
}
//+
Sphere(3) = {0, 1.5, 0, 2, -Pi/2, Pi/2, 2*Pi};
//+
Cylinder(4) = {-9.3, 1.5, 0.4, 10, 0, 0, 0.5, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/3} {
  Volume{4}; 
}
//+
BooleanUnion{ Volume{3}; Delete; }{ Volume{4}; Delete; }
//+
BooleanUnion{ Volume{3}; Delete; }{ Volume{1}; Delete; }
//+
BooleanUnion{ Volume{3}; Delete; }{ Volume{2}; Delete; }
