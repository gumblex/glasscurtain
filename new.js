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


function DrawLine(coord){
	if (coord[0]>=padding*1.5 && coord[0]<w-padding*1.5
		&& coord[1]>padding && coord[1]<=h - padding){
		//mouseDownOnElement=mouse;
		$("#tpoint").attr("transform","translate("+coord[0]+","+coord[1]+")");
		$("#tptext").attr("x", coord[0]+7);
		$("#tptext").attr("y", coord[1]+10);
		var nowDate = new Date();
		var startD = new Date(2014, 0);
		var nowDate = new Date(startD.setDate(Math.floor(yrevScale(coord[1]))+1));
		nowDate.setMinutes(xrevScale(coord[0]));
		$("#tptext").text((nowDate.getMonth()+1) + "-" + nowDate.getDate() + " " + ("0" + nowDate.getHours()).slice(-2) + ":" + ("0" + nowDate.getMinutes()).slice(-2));
		if (typeof map === "undefined") {return};
		cla = markerc.getPosition().lat() * pi180;
		cln = markerc.getPosition().lng() * pi180;
		dla = markerd.getPosition().lat() * pi180;
		dln = markerd.getPosition().lng() * pi180;
		for (var x=0;x<markeral.length;x++){
			ala = markeral[x].getPosition().lat() * pi180;
			aln = markeral[x].getPosition().lng() * pi180;
			bla = markerbl[x].getPosition().lat() * pi180;
			bln = markerbl[x].getPosition().lng() * pi180;
			wallh = gwallh[x];
			val=Calc(Math.floor(yrevScale(coord[1])),xrevScale(coord[0])/60,ala,aln,bla,bln,cla,cln,dla,dln,wallh*3,$("#planeh").spinner("value")*3,$("#window").slider("value"),x);
			if (val>-1){
			pab1l[x].setOptions({fillColor: Palette[val], visible: true})
			}else{pab1l[x].setOptions({visible: false})}
		}
	}
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
	var pE = object[0];
	var pF = object[1];
	var H = object[2];
	var pG = {
		x: (pA.y - pE.y - (pA.x * (pB.y - pA.y) / (pB.x - pA.x)) + Math.tan(-alpha) * pE.x) / (Math.tan(-alpha) - (pB.y - pA.y) / (pB.x - pA.x)),
		y: Math.tan(-alpha) * (pA.y - pE.y - (pA.x * (pB.y - pA.y) / (pB.x - pA.x)) + Math.tan(-alpha) * pE.x) / (Math.tan(-alpha) - (pB.y - pA.y) / (pB.x - pA.x)) + pE.y};
	var pH = {
		x: (pA.y - pF.y - (pA.x * (pB.y - pA.y) / (pB.x - pA.x)) + Math.tan(-alpha) * pF.x) / (Math.tan(-alpha) - (pB.y - pA.y) / (pB.x - pA.x)),
		y: (Math.tan(-alpha) * (pA.y - pF.y - (pA.x * (pB.y - pA.y) / (pB.x - pA.x)) + Math.tan(-alpha) * pF.x)) / (Math.tan(-alpha) - (pB.y - pA.y) / (pB.x - pA.x)) + pF.y};
	if ((pG.x < pA.x || pG.x > pB.x) && (pH.x < pA.x || pH.x > pB.x))
		{return [];}

	var d1 = Math.sqrt(Math.pow(pE.x - pG.x, 2) + Math.pow(pE.y - pG.y, 2));
	var d2 = Math.sqrt(Math.pow(pF.x - pH.x, 2) + Math.pow(pF.y - pH.y, 2));
	var L1 = H / Math.tan(h0);
	var Hmax1 = Math.tan(h0) * (L1 - d1);
	var Hmax2 = Math.tan(h0) * (L1 - d2);
	if (Hmax1 < 0 && Hmax2 < 0) {return [];}
	var L2 = Hmax1 / Math.tan(h0);
	var L3 = Hmax2 / Math.tan(h0);
	pG.la = pE.la - (d1 * Math.cos(Math.PI - alpha)) * 2 / dearth;
	pG.ln = pE.ln + (d1 * Math.sin(Math.PI - alpha)) * 2 / dearth;
	pH.la = pF.la - (d2 * Math.cos(Math.PI - alpha)) * 2 / dearth;
	pH.ln = pF.ln + (d2 * Math.sin(Math.PI - alpha)) * 2 / dearth;
	var pI = {la: pG.la - (L2 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pG.ln + (L2 * Math.sin(Math.PI - mu)) * 2 / dearth};
	var pJ = {la: pH.la - (L3 * Math.cos(Math.PI - mu)) * 2 / dearth,
	          ln: pH.ln + (L3 * Math.sin(Math.PI - mu)) * 2 / dearth};
	// [[pA.la,pA.ln],[pA2.la,pA2.ln],[pB2.la,pB2.ln],[pB.la,pB.ln]]
	if (Hmax1 > 0 && Hmax2 > 0) {
		return [{x:pG.la,y:pG.ln},{x:pH.la,y:pH.ln},{x:pJ.la,y:pJ.ln},{x:pI.la,y:pI.ln}];
	} else if (Hmax1 < 0) {
		return [{x:pH.la,y:pH.ln},{x:pJ.la,y:pJ.ln},{x:pI.la,y:pI.ln}];
	} else {
		return [{x:pG.la,y:pG.ln},{x:pJ.la,y:pJ.ln},{x:pI.la,y:pI.ln}];
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
	return [{x:pK.la,y:pK.ln},{x:pL.la,y:pL.ln},{x:pN.la,y:pN.ln},{x:pM.la,y:pM.ln}];
}
function mergeandclip(subject, clippers) {
	var polygons = [subject], newpoly;
	var i = 0, j = 0;
	for (i=0;i<clippers.length;i++) {
		if (clippers[i].length === 0) {continue}
		newpoly = [];
		for (j=0;j<polygons.length;j++) {
			if (polygons[j].length === 0) {continue}
			newpoly.concat(clip_polygon(polygons[j], clippers[i], "difference"));
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
	var pOla = markeral[0].getPosition().lat() * pi180;
	var pOln = markeral[0].getPosition().lng() * pi180;
	var pO = {la: pOla, ln: pOln};
	var Hp = $("#planeh").spinner("value")*3;
	var lw = $("#window").slider("value");
	var tla = 0, tln = 0;
	var avgla = 0, avgln = 0;
	var n = 0;
	for (i=0;i<markeral.length;i++){
		tla = markeral[i].getPosition().lat() * pi180;
		tln = markeral[i].getPosition().lng() * pi180;
		avgla += tla; avgln += tln;
		pA = {la: tla, ln: tln, x: recx(pOla,0,tla,0), y: recy(0,pOln,tla,tln)};
		if (walllength(pO, pA) > 1500) {return null;}
		tla = markerbl[i].getPosition().lat() * pi180;
		tln = markerbl[i].getPosition().lng() * pi180;
		avgla += tla; avgln += tln;
		pB = {la: tla, ln: tln, x: recx(pOla,0,tla,0), y: recy(0,pOln,tla,tln)};
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
	var pC = {la: tla, ln: tln, x: recx(pOla,0,tla,0), y: recy(0,pOln,tla,tln)};
	if (walllength(pO, pC) > 1500) {return null;}
	tla = markerd.getPosition().lat() * pi180;
	tln = markerd.getPosition().lng() * pi180;
	var pD = {la: tla, ln: tln, x: recx(pOla,0,tla,0), y: recy(0,pOln,tla,tln)};
	if (walllength(pO, pD) > 1500) {return null;}
	if (walllength(pC, pD) > 75) {return null;}
	var m = 0;
	if (pC.la !== pD.la) {
		m = Math.acos(Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)),2))/2)/Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)*Math.cos(Math.abs(pC.ln-pD.ln))),2)+Math.pow((Math.cos(pC.la)*Math.sin(Math.abs(pC.ln-pD.ln))),2))/2));
	}
	if (pC.la < pD.la) {m = Math.PI-m}
	if (pC.ln < pD.ln) {m = -m}
	var pE = {x: (pD.x*lw+pC.x*(100-lw))/100, y: (pD.y*lw+pC.y*(100-lw))/100};
	return [pAvg, walls, pC, pD, pE, m, Hp];
}
function lineGraph(index, range){
	if (visiblemark[index]){
		var i = 0, j = 0;
		var paths = [], newpath = [];
		for (i=0;i<range.length;i++) {
			newpath = [];
			for (j=0;j<range[i].length;j++) {
				newpath.push(new google.maps.LatLng(range[i][j].x/pi180, range[i][j].y/pi180));
			}
			paths.push(newpath);
		}
		pab1l[index].setOptions({
			paths: paths,
			visible: false,
		});
	}
}

function Calc(daynum,time,pAvg,walls,pC,pD,pE,m,Hp,draw) {
	var localtime = time+(pAvg.ln/pi180-120)/15;
	var a = Math.PI*2*(daynum+localtime/24)/365;
	var phi = 0.006918-0.399912*Math.cos(a)+0.070257*Math.sin(a)-0.006758*Math.cos(2*a)+0.000907*Math.sin(2*a)-0.002697*Math.cos(3*a)+0.001480*Math.sin(3*a);
	var h0=Math.asin(Math.cos((localtime-12)*15*pi180)*Math.cos(pAvg.la)*Math.cos(phi)+Math.sin(pAvg.la)*Math.sin(phi));
	if (h0 <= 0) {return -1;}
	var cosA = (Math.sin(h0)*Math.sin(pAvg.la)-Math.sin(phi))/Math.cos(h0)/Math.cos(pAvg.la);
	var alpha = Math.acos(cosA);
	if (localtime < 12) {alpha = -alpha};
	var pA, pB, n, L, lightrange;
	var worst = 0;
	for (var i=0;i<walls.length;i++) {
		pA = walls[i][0];
		pB = walls[i][1];
		L = walls[i][2]/Math.tan(h0);
		n = walls[i][3];
		// 阳光照在背面
		// if (Math.sin(alpha-n)>0) {if(draw !== -1) {return -2} else {return 0};}
		if (Math.sin(alpha-n)>0) {continue}
		var mu = 2*n-alpha;
		lightrange = doclips(h0, i, walls, [pC, pD, Hp], alpha, mu);
		if (draw !== -1) {lineGraph(draw, lightrange);}
		// 反射光照在背面
		if (Math.sin(mu-m)>0) {continue}
		if (pA.ln>pB.ln) {pB.y=-pB.y}; if (pA.la<pB.la) {pB.x=-pB.x};
		if (pA.ln>pC.ln) {pC.y=-pC.y}; if (pA.la<pC.la) {pC.x=-pC.x};
		if (pA.ln>pD.ln) {pD.y=-pD.y}; if (pA.la<pD.la) {pD.x=-pD.x};
		pF = {x: (pE.y-Math.tan(Math.PI-mu))/(pB.y/pB.x-Math.tan(Math.PI-mu)),
			  y: (pE.y-Math.tan(Math.PI-mu))*pB.y/((pB.y/pB.x-Math.tan(Math.PI-mu))*pB.x)};
		d = Math.sqrt(Math.pow(pE.x-pF.x,2)+Math.pow(pE.y-pF.y,2));
		hmax = Math.tan(h0)*(L-d);
		// 照不到高度
		if (hmax<Hp) {continue};
		xa=L*Math.cos(Math.PI-mu);
		xb=pB.x+L*Math.cos(Math.PI-mu);
		ya=L*Math.sin(Math.PI-mu);
		yb=pB.y+L*Math.sin(Math.PI-mu);
		if (!pointinrange(pE, lightrange)) {continue};

		gamma=-Math.PI/2+m-2*n+alpha;
		lambda=Math.acos(Math.abs(Math.cos(gamma)*Math.cos(h0)));
		B=rho*137000*Math.exp(-0.223/Math.sin(h0))/Math.PI;

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
	if (worst!==0) {console.log(worst,daynum,time)}
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

var pAvg = {la:0.5447030145459243,ln:2.1195261780590657};
var walls = [[{la:0.5447057772095524,ln:2.119526833446365,x:0,y:0},{la:0.5447002518822962,ln:2.119525522671766,x:35.21178896124298,y:7.144432904604358},60,0.20018106791599732]];
var pC = {la:0.5447036952043522,ln:2.119534323586931,x:13.26819649373995,y:40.82524573060971};
var pD = {la:0.5447000917276071,ln:2.11953404270666,x:36.23242228384485,y:39.29438480125423}
var pE = {x:24.7503093887924,y:40.059815265931974};
var m = 0.06656816376053089;
var Hp = 18;
var daynum = 123, time = 10;
for (daynum=0;daynum<365;daynum+=5) {
	for (time=0;time<1440;time+=6) {
		var localtime = time+(pAvg.ln/pi180-120)/15;
		var a = Math.PI*2*(daynum+localtime/24)/365;
		var phi = 0.006918-0.399912*Math.cos(a)+0.070257*Math.sin(a)-0.006758*Math.cos(2*a)+0.000907*Math.sin(2*a)-0.002697*Math.cos(3*a)+0.001480*Math.sin(3*a);
		var h0 = Math.asin(Math.cos((localtime-12)*15*pi180)*Math.cos(pAvg.la)*Math.cos(phi)+Math.sin(pAvg.la)*Math.sin(phi));
		var cosA = (Math.sin(h0)*Math.sin(pAvg.la)-Math.sin(phi))/Math.cos(h0)/Math.cos(pAvg.la);
		var alpha = Math.acos(cosA);
		if (localtime < 12) {alpha = -alpha};
		var n = 0.20018106791599732;
		var pA, pB, n, L, lightrange;
		var worst = 0;
		// 阳光照在背面
		// if (Math.sin(alpha-n)>0) {if(draw !== -1) {return -2} else {return 0};}
		if (Math.sin(alpha-n)>0) {print('continue')} else {
		var mu = 2*n-alpha;
		lightrange = doclips(h0, 0, walls, [pC, pD, Hp], alpha, mu);
		//print(JSON.stringify(clipshadow(h0, walls[0], [pC, pD, Hp], alpha, mu)));
		//print(JSON.stringify(lightrange));
}}}
