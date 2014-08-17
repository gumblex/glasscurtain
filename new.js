/*        /B
 *       /
 * 
 *    |A       \
 *    |         \
 *    |   \      \受照
 *         \
 *          \C
 * 
 * r =   入射光遮挡   玻璃幕墙, 遮挡墙面
 * 
 * r = clipshadow(h0, <墙A>, <受照>) + 
 *     clipshadow(h0, <墙A>, <墙B>) +
 *     clipshadow(h0, <墙A>, <墙C>)
 * cliplight(r, <墙B>) + cliplight(r, <墙C>)  ---.
 *                                              | 
 * r = clipshadow(h0, <墙B>, <受照>) +           |
 *     clipshadow(h0, <墙B>, <墙A>) +            |
 *     clipshadow(h0, <墙B>, <墙C>)              |
 * cliplight(r, <墙A>) + cliplight(r, <墙C>)  ---+>  --> 合并 & return？
 *                                              |      （得到总反射投影）
 * r = clipshadow(h0, <墙C>, <受照>) +           |
 *     clipshadow(h0, <墙C>, <墙A>) +            |
 *     clipshadow(h0, <墙C>, <墙B>)              |
 * cliplight(r, <墙A>) + cliplight(r, <墙B>)  ---'
 *                      
 * 
**/
'use strict';
var Palette=["#9edae5","#AEFF5C","#FFD34D","#ff7f0e","#e377c2","#d62728"];
var efname=["无","可接受","轻微影响","有影响","强影响","严重影响"];
var dearth = 12745594;
var pi180 = 0.0174532925199;
var timeoutID = 0;
var h = 320, w=400;
var padding = 20;
var lineHeight = 30;
var rho = 0.15;
var midleTime = 0;
var idleInterval = -1;
var mkdelTimeout = -1;
var xScale, yScale, xrevScale, yrevScale;
var nowcoord=[0,0];
var mouseDownOnElement=false;

// 经纬转直角
function recx(la0,_ln0,lax,_lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.sin(lax)-Math.sin(la0),2)+Math.pow(Math.cos(lax)-Math.cos(la0),2))/2)};
function recy(_la0,ln0,lax,lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.cos(lax)-Math.cos(lax)*Math.cos(Math.abs(ln0-lnx)),2)+Math.pow(Math.cos(lax)*Math.sin(Math.abs(ln0-lnx)),2))/2)};

function Vertex(vertex) {
	// Node in a circular doubly linked list.
	var self = this;
	self.la = vertex.la;
	self.ln = vertex.ln;
	self.next = null;
	self.prev = null;
	self.neighbour = null;
	self.entry = true;
	self.alpha = 0.0;
	self.intersect = false;
	self.checked = false;
}
function Vertex(vertex) {
	// Node in a circular doubly linked list.
	var self = this;
	self.la = vertex.la;
	self.ln = vertex.ln;
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
	var infinity = new Vertex({la: 1000000, ln: self.ln});
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
		p.push({la: v.la, ln: v.ln});
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
	var den = (c2.ln - c1.ln) * (s2.la - s1.la) - (c2.la - c1.la) * (s2.ln - s1.ln);
	if (!den) {return;}
	var us = ((c2.la - c1.la) * (s1.ln - c1.ln) - (c2.ln - c1.ln) * (s1.la - c1.la)) / den;
	var uc = ((s2.la - s1.la) * (s1.ln - c1.ln) - (s2.ln - s1.ln) * (s1.la - c1.la)) / den;
	if (((us === 0 || us === 1) && (0 <= uc && uc <= 1)) ||
	   ((uc === 0 || uc === 1) && (0 <= us && us <= 1))) {
		return;
	} else if ((0 < us && us < 1) && (0 < uc && uc < 1)) {
		var la = s1.la + us * (s2.la - s1.la);
		var ln = s1.ln + us * (s2.ln - s1.ln);
		return [{la: la, ln: ln}, us, uc];
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

// 检查顺时针，现在不需要
function checkclockwise(points) {
	var sum = 0;
	for (var i=1;i<points.length;i++) {
		sum += (points[i][0]-points[i-1][0])*(points[i][1]+points[i-1][1]);}
	sum += (points[i-1][0]-points[0][0])*(points[i-1][1]+points[0][1])
	if (sum < 0) {return points.reverse()} else {return points};
}

function zip(arrays) {
	return arrays[0].map(function(_,i){
		return arrays.map(function(array){return array[i]})
	});
}
function walllength(pA, pB) {
	var k = Math.sqrt(Math.pow(Math.sin(pB.la) - Math.sin(pA.la), 2) + Math.pow(Math.cos(pB.la) - Math.cos(pA.la) * Math.cos(Math.abs(pA.ln - pB.ln)), 2) + Math.pow(Math.cos(pA.la) * Math.sin(Math.abs(pA.ln - pB.ln)), 2));
	return dearth * Math.asin(k / 2);
}


// y = Math.tan(-alpha) * (x - x4) + y4;
// y = Math.tan(-alpha) * (x - x5) + y5;
// y = (y7 - y6) / (x7 - x6) * (x - x6) + y6;

// 入射光遮挡
// , pA2, pB2
// function clipshadow(h0, pA, pB, pE, pF, H) {
function clipshadow(h0, subject, object, alpha, mu) {
	// require {la, ln, x, y}
	var pA = subject[0];
	var pB = subject[1];
	var n = subject[2];
	var pE = object[0];
	var pF = object[1];
	var H = object[2];
	if (pA.la === pE.la) {
		var thetaea = Math.PI/2;
		if (pA.ln < pE.ln) {thetaea = -thetaea} else {return []}
	} else {
		var thetaea = Math.atan((pA.ln-pE.ln)/(pA.la-pE.la));
		if (pA.la < pE.la) {thetaea = Math.PI - thetaea}
	}
	if (pA.la === pF.la) {
		var thetafa = Math.PI/2;
		if (pA.ln < pF.ln) {thetafa = -thetafa} else {return []}
	} else {
		var thetafa = Math.atan((pA.ln-pF.ln)/(pA.la-pF.la));
		if (pA.la < pF.la) {thetaea = Math.PI - thetafa}
	}
	if (Math.sin(thetaea-n) > 0 && Math.sin(thetafa-n) > 0) {return []}
	var pG = {
		la: (pA.ln - pE.ln - (pA.la * (pB.ln - pA.ln) / (pB.la - pA.la)) + Math.tan(-alpha) * pE.la) / (Math.tan(-alpha) - (pB.ln - pA.ln) / (pB.la - pA.la)),
		ln: Math.tan(-alpha) * ((pA.ln - pE.ln - (pA.la * (pB.ln - pA.ln) / (pB.la - pA.la)) + Math.tan(-alpha) * pE.la) / (Math.tan(-alpha) - (pB.ln - pA.ln) / (pB.la - pA.la)) - pE.la) + pE.ln};
	var pH = {
		la: (pA.ln - pF.ln - (pA.la * (pB.ln - pA.ln) / (pB.la - pA.la)) + Math.tan(-alpha) * pF.la) / (Math.tan(-alpha) - (pB.ln - pA.ln) / (pB.la - pA.la)),
		ln: Math.tan(-alpha) * ((pA.ln - pF.ln - (pA.la * (pB.ln - pA.ln) / (pB.la - pA.la)) + Math.tan(-alpha) * pF.la) / (Math.tan(-alpha) - (pB.ln - pA.ln) / (pB.la - pA.la)) - pF.la) + pF.ln};
	if ((pG.la < pA.la || pG.la > pB.la) && (pH.la < pA.la || pH.la > pB.la))
		{return []}

	//var d1 = Math.sqrt(Math.pow(pE.x - pG.x, 2) + Math.pow(pE.y - pG.y, 2));
	//var d2 = Math.sqrt(Math.pow(pF.x - pH.x, 2) + Math.pow(pF.y - pH.y, 2));
	var d1 = walllength(pE,pG);
	var d2 = walllength(pF,pH);
	var L1 = H / Math.tan(h0);
	var Hmax1 = Math.tan(h0) * (L1 - d1);
	var Hmax2 = Math.tan(h0) * (L1 - d2);
	if (Hmax1 < 0 && Hmax2 < 0) {return []}
	var L2 = L1 - d1;
	var L3 = L1 - d2;
	//pG.la = pE.la - (d1 * Math.cos(Math.PI - alpha)) * 2 / dearth;
	//pG.ln = pE.ln + (d1 * Math.sin(Math.PI - alpha)) * 2 / dearth;
	//pH.la = pF.la - (d2 * Math.cos(Math.PI - alpha)) * 2 / dearth;
	//pH.ln = pF.ln + (d2 * Math.sin(Math.PI - alpha)) * 2 / dearth;
	var pI = {la: pG.la - (L2 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pG.ln + (L2 * Math.sin(Math.PI - mu)) * 2 / dearth};
	var pJ = {la: pH.la - (L3 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pH.ln + (L3 * Math.sin(Math.PI - mu)) * 2 / dearth};
	// [[pA.la,pA.ln],[pA2.la,pA2.ln],[pB2.la,pB2.ln],[pB.la,pB.ln]]
	if (Hmax1 > 0 && Hmax2 > 0) {
		return [pG,pH,pJ,pI];
	} else if (Hmax1 < 0) {
		return [pH,pJ,pI];
	} else {
		return [pG,pJ,pI];
	}
}
function cliplight(h0, object, mu) {
	var pK = object[0];
	var pL = object[1];
	var H = object[2];
	var L4 = H / Math.tan(h0);
	var pM = {la: pK.la - (L4 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pK.ln + (L4 * Math.sin(Math.PI - mu)) * 2 / dearth};
	var pN = {la: pL.la - (L4 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pL.ln + (L4 * Math.sin(Math.PI - mu)) * 2 / dearth};
	return [pK,pL,pN,pM];
}
function mergeandclip(subject, clippers) {
	var polygons = [subject], newpoly;
	var i = 0, j = 0;
	for (i=0;i<clippers.length;i++) {
		if (clippers[i].length === 0) {continue}
		newpoly = [];
		for (j=0;j<polygons.length;j++) {
			if (polygons[j].length !== 0) {
			newpoly = newpoly.concat(clip_polygon(polygons[j], clippers[i], "difference"));
			}
		}
		polygons = newpoly;
	}
	return polygons;
}
function doclips(h0, gindex, walls, object, alpha, mu) {
	var j = 0;
	var cs = clipshadow(h0, walls[gindex], object, alpha, mu);
	var cl = [];
	if (cs.length !== 0) {cl.push(cs)}
	for (j=0;j<walls.length;j++) {
		if (j !== gindex) {
			cs = clipshadow(h0, walls[gindex], walls[j], alpha, mu);
			if (cs.length !== 0) {cl.push(cs)}
			cl.push(cliplight(h0, walls[j], mu))
		}
	}
	var origpoly = [walls[gindex][0],
	 {la: walls[gindex][0].la-(walls[gindex][2]*Math.cos(Math.PI-mu)*2/dearth),
	  ln: walls[gindex][0].ln+(walls[gindex][2]*Math.sin(Math.PI-mu)*2/dearth)},
	 {la: walls[gindex][1].la-(walls[gindex][2]*Math.cos(Math.PI-mu)*2/dearth),
	  ln: walls[gindex][1].ln+(walls[gindex][2]*Math.sin(Math.PI-mu)*2/dearth)},
	 walls[gindex][1]];
	if (cl.length === 0) {return [origpoly]}
	else {return mergeandclip(origpoly, cl)}
}
function pointInPolygon(la,ln,poly) {
	var polySides = poly.length;
	var i = 0, j = polySides-1;
	var oddNodes = false;
	for (i=0; i<polySides; i++) {
		if ((poly[i].ln< ln && poly[j].ln>=ln
		   ||poly[j].ln< ln && poly[i].ln>=ln)
		   &&(poly[i].la<=la || poly[j].la<=la)) {
			if (poly[i].la+(ln-poly[i].ln)/(poly[j].ln-poly[i].ln)*(poly[j].la-poly[i].la)<la) {oddNodes=!oddNodes}
		}
		j = i;
	}
	return oddNodes;
}

function pointinrange(point, range) {
	for (var i=0;i<range.length;i++) {
		// xyarray = zip(range[i])
		if (pointInPolygon(point.x,point.y,range[i])) {
			return true;}
	}
	return false;
}

function PrepareData() {
	var i = 0;
	var walls = [];
	var pA, pB;
	var pO = {la: markeral[0].getPosition().lat() * pi180,
	          ln: markeral[0].getPosition().lng() * pi180};
	var Hp = $("#planeh").spinner("value")*3;
	var lw = $("#window").slider("value");
	var tla = 0, tln = 0;
	var avgla = 0, avgln = 0;
	var n = 0;
	for (i=0;i<markeral.length;i++){
		tla = markeral[i].getPosition().lat() * pi180;
		tln = markeral[i].getPosition().lng() * pi180;
		avgla += tla; avgln += tln;
		pA = {la: tla, ln: tln,
		      x: recx(pO.la,0,tla,0), y: recy(0,pO.ln,tla,tln)};
		if (pO.ln>pA.ln) {pA.y=-pA.y}; if (pO.la<pA.la) {pA.x=-pA.x};
		if (walllength(pO, pA) > 1500) {return null;}
		tla = markerbl[i].getPosition().lat() * pi180;
		tln = markerbl[i].getPosition().lng() * pi180;
		avgla += tla; avgln += tln;
		pB = {la: tla, ln: tln,
		      x: recx(pO.la,0,tla,0), y: recy(0,pO.ln,tla,tln)};
		if (pO.ln>pB.ln) {pB.y=-pB.y}; if (pO.la<pB.la) {pB.x=-pB.x};
		if (walllength(pO, pB) > 1500) {return null;}
		if (walllength(pA, pB) > 200) {return null;}
		n = 0;
		if (pA.la !== pB.la) {
			n = Math.acos(Math.asin(Math.sqrt(Math.pow((Math.sin(pB.la)-Math.sin(pA.la)),2)+Math.pow((Math.cos(pB.la)-Math.cos(pA.la)),2))/2)/Math.asin(Math.sqrt(Math.pow((Math.sin(pB.la)-Math.sin(pA.la)),2)+Math.pow((Math.cos(pB.la)-Math.cos(pA.la)*Math.cos(Math.abs(pA.ln-pB.ln))),2)+Math.pow((Math.cos(pA.la)*Math.sin(Math.abs(pA.ln-pB.ln))),2))/2));
		}
		if (pA.la < pB.la) {n = Math.PI-n}
		if (pA.ln < pB.ln) {n = -n}
		walls.push([pA, pB, gwallh[i]*3, n]);
	}
	var pAvg = {la: avgla/markeral.length/2, ln: avgln/markeral.length/2};
	tla = markerc.getPosition().lat() * pi180;
	tln = markerc.getPosition().lng() * pi180;
	var pC = {la: tla, ln: tln,
	          x: recx(pO.la,0,tla,0), y: recy(0,pO.ln,tla,tln)};
	if (pO.ln>pC.ln) {pC.y=-pC.y}; if (pO.la<pC.la) {pC.x=-pC.x};
	if (walllength(pO, pC) > 1500) {return null;}
	tla = markerd.getPosition().lat() * pi180;
	tln = markerd.getPosition().lng() * pi180;
	var pD = {la: tla, ln: tln,
	          x: recx(pO.la,0,tla,0), y: recy(0,pO.ln,tla,tln)};
	if (pO.ln>pD.ln) {pD.y=-pD.y}; if (pO.la<pD.la) {pD.x=-pD.x};
	if (walllength(pO, pD) > 1500) {return null;}
	if (walllength(pC, pD) > 75) {return null;}
	var m = 0;
	if (pC.la !== pD.la) {
		m = Math.acos(Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)),2))/2)/Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)*Math.cos(Math.abs(pC.ln-pD.ln))),2)+Math.pow((Math.cos(pC.la)*Math.sin(Math.abs(pC.ln-pD.ln))),2))/2));
	}
	if (pC.la < pD.la) {m = Math.PI-m}
	if (pC.ln < pD.ln) {m = -m}
	var pE = {la: (pD.la*lw+pC.la*(100-lw))/100,
	          ln: (pD.ln*lw+pC.ln*(100-lw))/100};
	return [pAvg, walls, pC, pD, pE, pO, m, Hp];
}
function lineGraph(pO, index, range){
	if (visiblemark[index] && !pab1l[index].get('invisible')){
		var i = 0, j = 0;
		var paths = [], newpath = [];
		for (i=0;i<range.length;i++) {
			newpath = [];
			for (j=0;j<range[i].length;j++) {
				if (!range[i][j].hasOwnProperty('la')) {
					range[i][j].la = pO.la - range[i][j].x*2/dearth;
					range[i][j].ln = pO.ln + range[i][j].y*2/dearth;
				}
				newpath.push(new google.maps.LatLng(range[i][j].la/pi180, range[i][j].ln/pi180));
			}
			paths.push(newpath);
		}
		pab1l[index].setOptions({
			paths: paths,
			visible: false,
		});
	}
}

function Calc(daynum,time,pAvg,walls,pC,pD,pE,pO,m,Hp,draw) {
	var localtime = time+(pAvg.ln/pi180-120)/15;
	var a = Math.PI*2*(daynum+localtime/24)/365;
	var phi = 0.006918-0.399912*Math.cos(a)+0.070257*Math.sin(a)-0.006758*Math.cos(2*a)+0.000907*Math.sin(2*a)-0.002697*Math.cos(3*a)+0.001480*Math.sin(3*a);
	var h0 = Math.asin(Math.cos((localtime-12)*15*pi180)*Math.cos(pAvg.la)*Math.cos(phi)+Math.sin(pAvg.la)*Math.sin(phi));
	if (h0 <= 0) {return -1}
	var cosA = (Math.sin(h0)*Math.sin(pAvg.la)-Math.sin(phi))/Math.cos(h0)/Math.cos(pAvg.la);
	var alpha = Math.acos(cosA);
	if (localtime < 12) {alpha = -alpha};
	var pA, pB, pF, n, L, mu, lightrange, d, hmax, gamma, lambda, B;
	var worst = 0;
	for (var i=0;i<walls.length;i++) {
		pA = walls[i][0];
		pB = walls[i][1];
		L = walls[i][2]/Math.tan(h0);
		n = walls[i][3];
		// 阳光照在背面
		// if (Math.sin(alpha-n)>0) {if(draw !== -1) {return -2} else {return 0};}
		
		if (Math.sin(alpha-n)>0) {
			if (draw) {pab1l[i].set('invisible', true)}
			continue;
		}
		mu = 2*n-alpha;
		if (!draw) {if (Math.sin(mu-m)>0) {continue}}
		// 反射光照在背面
		lightrange = doclips(h0, i, walls, [pC, pD, Hp], alpha, mu);
		// console.log(JSON.stringify(lightrange));
		if (draw) {
			pab1l[i].set('invisible', false);
			lineGraph(pO, i, lightrange);
			if (Math.sin(mu-m)>0) {continue}
		}
		if (pA.ln>pD.ln) {pD.y=-pD.y}; if (pA.la<pD.la) {pD.x=-pD.x};
		pF = {la: (pE.ln-Math.tan(Math.PI-mu))/(pB.ln/pB.la-Math.tan(Math.PI-mu)),
		      ln: (pE.ln-Math.tan(Math.PI-mu))*pB.ln/((pB.ln/pB.la-Math.tan(Math.PI-mu))*pB.la)};
		d = walllength(pE, pF);
		hmax = Math.tan(h0)*(L-d);
		// 照不到高度
		if (hmax<Hp) {continue};
		if (!pointinrange(pE, lightrange)) {continue};

		gamma = -Math.PI/2+m-2*n+alpha;
		lambda = Math.acos(Math.abs(Math.cos(gamma)*Math.cos(h0)));
		B = rho*137000*Math.exp(-0.223/Math.sin(h0))/Math.PI;

		if (lambda<Math.PI/12) {
			if (B>=2000) {if (worst<5) {worst = 5}; continue;};
		} else if (lambda>=Math.PI/12 && lambda<=Math.PI/6) {
			if (B>=2000 && B<4000) {
				if (worst<2) {worst = 2}; continue;
			} else if (B>=4000 && B<6000) {
				if (worst<3) {worst = 3}; continue;
			} else if (B>=6000) {
				if (worst<4) {worst = 4}; continue;
			};
		} else if (lambda>Math.PI/6) {
			if (worst<1) {worst = 1}; continue;
		};
		// continue;
	}
	return worst;
}

/*
var i = 0, j = 0, val = 0;
for (i=0;i<365;i+=5) {
	for (j=0;j<1440;j+=6) {
		val = Calc(i,j/60,{la:0.5447030145459243,ln:2.1195261780590657},[[{la:0.5447057772095524,ln:2.119526833446365,x:0,y:0},{la:0.5447002518822962,ln:2.119525522671766,x:35.21178896124298,y:7.144432904604358},60,0.20018106791599732]],{la:0.5447036952043522,ln:2.119534323586931,x:13.26819649373995,y:40.82524573060971},{la:0.5447000917276071,ln:2.11953404270666,x:36.23242228384485,y:39.29438480125423},{x:24.7503093887924,y:40.059815265931974},0.06656816376053089,18,-1)
		//if (val===0) {
		//print(val+','+i+','+j)}
	}
}*/

var CalcPara = [{la:0.5451561292589383,ln:2.1196476880399873},[[{la:0.5451519063426971,ln:2.119644569098642,x:0,y:0},{la:0.5451573100747598,ln:2.1196475651548687,x:-34.43688747800468,y:16.32560879674119},60,-2.698900060621317],[{la:0.545157690336719,ln:2.1196501632973774,x:-36.86021975060163,y:30.482965643321073},{la:0.5451576102815774,ln:2.119648454609061,x:-36.35004458521584,y:21.172270052390783},60,1.5160565471943037]],{la:0.5451540278099882,ln:2.1196491099963604,x:-13.519680388541083,y:24.74355008465445},{la:0.5451553487222061,ln:2.119649344063253,x:-21.93758580766993,y:26.018969966558025},{la:0.5451546882660971,ln:2.1196492270298064},{la:0.5451519063426971,ln:2.119644569098642},-2.991221225676089,57];
var pAvg = CalcPara[0];
var walls = CalcPara[1];
var pC = CalcPara[2];
var pD = CalcPara[3];
var pE = CalcPara[4];
var pO = CalcPara[5];
var m = CalcPara[6];
var Hp = CalcPara[7];

var daynum = 123, time = 10;
/*for (daynum=0;daynum<365;daynum++) {
	for (time=0;time<1440;time+=2) {
		var localtime = time+(pAvg.ln/pi180-120)/15;
		var a = Math.PI*2*(daynum+localtime/24)/365;
		var phi = 0.006918-0.399912*Math.cos(a)+0.070257*Math.sin(a)-0.006758*Math.cos(2*a)+0.000907*Math.sin(2*a)-0.002697*Math.cos(3*a)+0.001480*Math.sin(3*a);
		var h0 = Math.asin(Math.cos((localtime-12)*15*pi180)*Math.cos(pAvg.la)*Math.cos(phi)+Math.sin(pAvg.la)*Math.sin(phi));
		var cosA = (Math.sin(h0)*Math.sin(pAvg.la)-Math.sin(phi))/Math.cos(h0)/Math.cos(pAvg.la);
		var alpha = Math.acos(cosA);
		if (localtime < 12) {alpha = -alpha};
		var n = 0.2547583024269751;
		var pA, pB, n, L, lightrange;
		var worst = 0;
		var val = Calc(daynum,time/60,pAvg,walls,pC,pD,pE,pO,m,Hp,false);
		if (val!==-1){print(daynum,time,val)};
}}*/
for (time=900;time<1000;time++) {
//print(JSON.stringify([pAvg,walls,pC,pD,pE,pO,m,Hp]));
print(Calc(126,time/60,pAvg,walls,pC,pD,pE,pO,m,Hp,false))
}
/*
110 328 1
110 330 0
110 332 1
110 334 0
110 336 1
110 338 0
110 340 1
*/
