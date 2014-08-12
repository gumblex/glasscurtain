/* Efficient Clipping of Arbitrary Polygons
Based on the paper "Efficient Clipping of Arbitrary Polygons" by Günther
Greiner (greiner[at]informatik.uni-erlangen.de) and Kai Hormann
(hormann[at]informatik.tu-clausthal.de), ACM Transactions on Graphics
1998;17(2):71-83.
Available at: http://www.inf.usi.ch/hormann/papers/Greiner.1998.ECO.pdf
You should have received the README file along with this program.
If not, see <https://github.com/helderco/polyclip>
*/
function Vertex(vertex) {
	// Node in a circular doubly linked list.
	var self = this;
	self.x = vertex.x;
	self.y = vertex.y;
	self.next = null;
	self.prev = null;
	self.neighbour = null;
	self.entry = true;
	self.alpha = 0.0;
	self.intersect = false;
	self.checked = false;
}
Vertex.prototype.isInside = function (poly) {
	// Test if a vertex lies inside a polygon (odd-even rule).
	var self = this;
	var winding_number = 0;
	var infinity = new Vertex({x: 1000000, y: self.y});
	var polyiter = poly.iter(), q, kq = 0;
	for (kq=0;kq<polyiter.length;kq++) {
		q = polyiter[kq];
		if (!q.intersect && !!intersect(self, infinity, q, poly.next(q.next))) {
			winding_number++;}
	}
	return ((winding_number % 2) !== 0);
};
Vertex.prototype.setChecked = function () {
	var self = this;
	self.checked = true;
	if (self.neighbour) {
		if (!self.neighbour.checked) {self.neighbour.setChecked();}
	}
};
function Polygon() {this.first = null;}
Polygon.prototype.add = function (vertex) {
	var self = this;
	if (!self.first) {
		self.first = vertex;
		self.first.next = vertex;
		self.first.prev = vertex;
	} else {
		var next = self.first;
		var prev = next.prev;
		next.prev = vertex;
		vertex.next = next;
		vertex.prev = prev;
		prev.next = vertex;
	}
};
Polygon.prototype.insert = function (vertex, start, end) {
	// Insert and sort a vertex between a specified pair of vertices.
	var curr = start;
	while (curr !== end && curr.alpha < vertex.alpha) {curr = curr.next;}
	vertex.next = curr;
	var prev = curr.prev;
	vertex.prev = prev;
	prev.next = vertex;
	curr.prev = vertex;
};
Polygon.prototype.next = function (v) {
	var c = v;
	while (c.intersect) {c = c.next;}
	return c;
};
Polygon.prototype.first_intersect = function () {
	var selfiter = this.iter(), v, kv = 0;
	for (kv=0;kv<selfiter.length;kv++) {
		v = selfiter[kv];
		if (v.intersect && !v.checked) {break;}
	}
	return v;
};
Polygon.prototype.points = function () {
	var p = [];
	var selfiter = this.iter(), v, kv = 0;
	for (kv=0;kv<selfiter.length;kv++) {
		v = selfiter[kv];
		p.push({x: v.x, y: v.y});
	}
	return p;
};
Polygon.prototype.unprocessed = function () {
	var selfiter = this.iter(), v, kv = 0;
	for (kv=0;kv<selfiter.length;kv++) {
		v = selfiter[kv];
		if (v.intersect && !v.checked) {return true;}
	}
	return false;
};
Polygon.prototype.clip = function (clip, s_entry, c_entry) {
	// Clip self polygon using another one as a clipper.
	var self = this;
	var selfiter = self.iter(), s, ks = 0;
	var clipiter = clip.iter(), c, kc = 0;
	var islist, iS, iC;
	for (ks=0;ks<selfiter.length;ks++) {
		s = selfiter[ks];
		if (!s.intersect) {
        clipiter = clip.iter();
        for (kc=0;kc<clipiter.length;kc++) {
            c = clipiter[kc];
            if (!c.intersect) {
                islist = intersect(s, self.next(s.next), c, clip.next(c.next));
                if (islist !== undefined) {
                    iS = new Vertex(islist[0]);
                    iS.alpha = islist[1];
                    iS.intersect = true;
                    iS.entry = false;
                    iC = new Vertex(islist[0]);
                    iC.alpha = islist[2];
                    iC.intersect = true;
                    iC.entry = false;
                    iS.neighbour = iC;
                    iC.neighbour = iS;
                    self.insert(iS, s, self.next(s.next));
                    clip.insert(iC, c, clip.next(c.next));
                }
            }
        }}
	}
	s_entry = (!s_entry !== !self.first.isInside(clip));
	selfiter = self.iter();
	for (ks=0;ks<selfiter.length;ks++) {
		s = selfiter[ks];
		if (s.intersect) {
			s.entry = s_entry;
			s_entry = !s_entry;
		}
	}
	c_entry = (!c_entry !== !clip.first.isInside(self));
	clipiter = clip.iter();
	for (kc=0;kc<clipiter.length;kc++) {
		c = clipiter[kc];
		if (c.intersect) {
			c.entry = c_entry;
			c_entry = !c_entry;
		}
	}
	var list = [], current, clipped;
	while (self.unprocessed()) {
		current = self.first_intersect();
		clipped = new Polygon();
		clipped.add(new Vertex(current));
		while (true) {
			current.setChecked();
			if (current.entry) {
				while (true) {
					current = current.next;
					clipped.add(new Vertex(current));
					if (current.intersect) {break;}
				}
			} else {
				while (true) {
					current = current.prev;
					clipped.add(new Vertex(current));
					if (current.intersect) {break;}
				}
			}
			current = current.neighbour;
			if (current.checked) {break;}
		}
		list.push(clipped);
	}
	if (list.length === 0) {list.push(self);}
	return list;
};
Polygon.prototype.iter = function () {
	var self = this;
	var s = self.first;
	var templist = [];
	while (true) {
		templist.push(s);
		s = s.next;
		if (s === self.first) {return templist;}
	}
};
function intersect(s1, s2, c1, c2) {
	var den = (c2.y - c1.y) * (s2.x - s1.x) - (c2.x - c1.x) * (s2.y - s1.y);
	if (!den) {return;}
	var us = ((c2.x - c1.x) * (s1.y - c1.y) - (c2.y - c1.y) * (s1.x - c1.x)) / den;
	var uc = ((s2.x - s1.x) * (s1.y - c1.y) - (s2.y - s1.y) * (s1.x - c1.x)) / den;
	if (((us === 0 || us === 1) && (0 <= uc && uc <= 1)) ||
	   ((uc === 0 || uc === 1) && (0 <= us && us <= 1))) {
		return;
	} else if ((0 < us && us < 1) && (0 < uc && uc < 1)) {
		var x = s1.x + us * (s2.x - s1.x);
		var y = s1.y + us * (s2.y - s1.y);
		return [{x: x, y: y}, us, uc];
	}
	return;
}
function clip_polygon(subject, clipper, operation) {
	var Subject = new Polygon(), Clipper = new Polygon();
	var k = 0;
	for (k=0;k<subject.length;k++) {
		Subject.add(new Vertex(subject[k]));}
	for (k=0;k<clipper.length;k++) {
		Clipper.add(new Vertex(clipper[k]));}
	var clipped;
	if (operation === 'difference') {
		clipped = Subject.clip(Clipper, false, true);
	} else if (operation === 'reversed-diff') {
		clipped = Clipper.clip(Subject, false, true);
	} else if (operation === 'union') {
		clipped = Subject.clip(Clipper, false, false);
	} else if (operation === 'intersection') {
		clipped = Subject.clip(Clipper, true, true);
	}
	var clippedlist = [];
	for (k=0;k<clipped.length;k++) {clippedlist.push(clipped[k].points());}
	return clippedlist;
}
