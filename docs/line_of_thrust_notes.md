# Line of Thrust Notes

OK, let's go back to the sweep_left, sweep_right idea. First of all, let's review some basics.

The Qi value computed by Spencer's method is based on the following assumptions:

Qi = sum of resultant side forces on the two sides of a slice. They act at an angle theta. They are all parallel. 

Also,

Qi = sum of forces on base of slice due to weight of slice, normal force, pore pressure, and mobilized shear resistance. We calculate Qi using these forces, but they still are assumed to be equal to sum of side forces on the two sides of the slice.

It would make sense that the sign of Qi might change as you move from the top of the slope to the bottom. At the top of the slope, the weight and the normal forces are pushing toward failure. But somewhere near the middle of the failure surface, the the weight and the normal forces push in the other direction and act to make the slope stable. So the sign switching is NOT violating any physical principles.

So to do the iterative sweep, let's say we are starting from the left and it is a right-facing slope. For the first slice, Q1 = some of both side forces. But the force on the left side is zero, so Q1 is applied to the right side. Once we have that, we can solve for the y_thrust value. Then we move to the second slice. We take side force we calculated on the right side of slice 1 (call this R_1,2) and apply it to the left side of slice 2. Then, to get the side force on the right side of slice 2 (R_2,3), we take Q2 for slice 2 and subtract the force we found from the first slice. In other words, R_2,3 = Q2 - R_1,2. We then proceed to the right. Theoretically, when we get all the way to the right, the forces should balance, but they might not exactly. So we could work to the center and do some kind of averaging.

Exactly right — and this is a principled, mechanics-based way to construct the line of thrust that handles the changing direction of interslice forces while staying within Spencer’s framework.

Let’s recap and structure your plan:

⸻

✅ Key Assumptions from Spencer’s Method
	•	Q_i is the net resultant of the two interslice side forces on slice i, acting at angle \theta.
	•	Spencer assumes:
	•	\vec{Z}{i-1,i} + \vec{Z}{i,i+1} = \vec{Q}_i
	•	All \vec{Z}’s are parallel (i.e., same direction θ)
	•	Q_i can change sign, and doing so is physically valid
	•	You’ve correctly noted:
	•	Side forces switch direction near the center of the slope due to change from driving to resisting forces.

⸻

🔁 Your Sweep Logic (Iterative Construction)

For a right-facing slope (starting from left):
	1.	Initialize:
	•	Left interslice force on slice 0 is zero
	•	So R_{1} = Q_0 acts entirely on right side of slice 0
	2.	For each slice i:
	•	Left interslice force: \vec{Z}_{i-1,i} = carried forward from previous slice
	•	Known \vec{Q}_i
	•	Compute right interslice force:
\vec{Z}_{i,i+1} = \vec{Q}i - \vec{Z}{i-1,i}
	3.	For each interface:
	•	Resolve the direction vector (angle θ)
	•	Moment balance gives vertical position (y) of each interslice force
	4.	At the end, compare net residual force. If needed, do:
	•	A second sweep from the right
	•	Average left/right values