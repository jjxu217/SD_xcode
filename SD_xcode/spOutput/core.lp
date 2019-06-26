\ENCODING=ISO-8859-1
\Problem name: ./spInput/SUFL/SUFL_F2_MinF1_FO1_D3_S4/SUFL_F2_MinF1_FO1_D3_S4.cor

Minimize
 OBJ: 30 x0 + 30 x1 + y0_0 + y1_0 + y2_0 + y0_1 + y1_1 + y2_1 + y0_2 + y1_2
      + y2_2
Subject To
 sum_of_facilities:             x0 + x1 >= 1
 sum_of_demand_0:               y0_0 + y1_0 + y2_0 >= 1
 sum_of_demand_1:               y0_1 + y1_1 + y2_1 >= 1
 sum_of_demand_2:               y0_2 + y1_2 + y2_2 >= 1
 facilities_demands_0_0:        - x0 + y0_0 <= 0
 facilities_demands_0_1:        - x0 + y0_1 <= 0
 facilities_demands_0_2:        - x0 + y0_2 <= 0
 facilities_demands_1_0:        - x1 + y1_0 <= 0
 facilities_demands_1_1:        - x1 + y1_1 <= 0
 facilities_demands_1_2:        - x1 + y1_2 <= 0
 opened_facilities_demands_2_0: y2_0 <= 1
 opened_facilities_demands_2_1: y2_1 <= 1
 opened_facilities_demands_2_2: y2_2 <= 1
Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= y0_0 <= 1
 0 <= y1_0 <= 1
 0 <= y2_0 <= 1
 0 <= y0_1 <= 1
 0 <= y1_1 <= 1
 0 <= y2_1 <= 1
 0 <= y0_2 <= 1
 0 <= y1_2 <= 1
 0 <= y2_2 <= 1
Binaries
 x0  x1  y0_0  y1_0  y2_0  y0_1  y1_1  y2_1  y0_2  y1_2  y2_2 
End
