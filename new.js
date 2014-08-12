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

// 经纬转直角
function recx(la0,_ln0,lax,_lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.sin(lax)-Math.sin(la0),2)+Math.pow(Math.cos(lax)-Math.cos(la0),2))/2)};
function recy(_la0,ln0,lax,lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.cos(lax)-Math.cos(lax)*Math.cos(Math.abs(ln0-lnx)),2)+Math.pow(Math.cos(lax)*Math.sin(Math.abs(ln0-lnx)),2))/2)};

// 检查顺时针，现在不需要
function checkclockwise(points) {
	var sum = 0;
	for (var i=1;i<points.length;i++) {
		sum += (points[i][0]-points[i-1][0])*(points[i][1]+points[i-1][1]);}
	sum += (points[i-1][0]-points[0][0])*(points[i-1][1]+points[0][1])
	if (sum < 0) {return points.reverse()} else {return points};
}

function allclips(h0, glasses, object) {
	var cl, ranges;
	var i = 0, j = 0;
	for (i=0;i<glasses.length;i++) {
		cl = [clipshadow(h0, glasses[i], object)]
		for (j=0;j<glasses.length;j++) {
			if (j !== i) {
				cl.push(clipshadow(h0, glasses[i], glasses[j]))
				cl.push(cliplight(glasses[j]))
			}
		}
		ranges.push(mergeandclip(glasses[i], cl))
	}
	return ranges
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
function clipshadow(h0, subject, object) {
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
	pG.la = pE.la - (d1 * Math.cos(Math.PI - alpha)) / r;
	pG.ln = pE.ln + (d1 * Math.sin(Math.PI - alpha)) / r;
	pH.la = pF.la - (d2 * Math.cos(Math.PI - alpha)) / r;
	pH.ln = pF.ln + (d2 * Math.sin(Math.PI - alpha)) / r;
	var pI = {la: pG.la - (L2 * Math.cos(Math.PI - mu)) / r,
	          ln: pG.ln + (L2 * Math.sin(Math.PI - mu)) / r};
	var pJ = {la: pH.la - (L3 * Math.cos(Math.PI - mu)) / r,
	          ln: pH.ln + (L3 * Math.sin(Math.PI - mu)) / r};
	// [[pA.la,pA.ln],[pA2.la,pA2.ln],[pB2.la,pB2.ln],[pB.la,pB.ln]]
	if (Hmax1 > 0 && Hmax2 > 0) {
		return [[pG.la,pG.ln],[pH.la,pH.ln],[pJ.la,pJ.ln],[pI.la,pI.ln]];
	} else if (Hmax1 < 0) {
		return [[pH.la,pH.ln],[pJ.la,pJ.ln],[pI.la,pI.ln]];
	} else {
		return [[pG.la,pG.ln],[pJ.la,pJ.ln],[pI.la,pI.ln]];
	}
}

// 反射光遮挡
function cliplight(h0, object, mu) {
	var pK = object[0];
	var pL = object[1];
	var H = object[2];
	var L4 = H / Math.tan(h0);
	var pM = {la: pK.la - (L4 * Math.cos(Math.PI - mu)) / r,
	          ln: pK.ln + (L4 * Math.sin(Math.PI - mu)) / r};
	var pN = {la: pL.la - (L4 * Math.cos(Math.PI - mu)) / r,
	          ln: pL.ln + (L4 * Math.sin(Math.PI - mu)) / r};
	return [[pK.la,pK.ln],[pL.la,pL.ln],[pN.la,pN.ln],[pM.la,pM.ln]];
}

function mergeandclip(subject, clippers) {
	var polygons = [subject], newpoly;
	var i = 0, j = 0;
	for (i=0;i<clippers.length;i++) {
		newpoly = [];
		for (j=0;j<polygons.length;j++) {
			newpoly.concat(clip_polygon(polygons[j], clippers[i], "difference"));
		}
		polygons = newpoly;
	}
	return polygons;
}

// 综合遮挡

function doclips(h0, gindex, glasses, object, mu) {
	var j = 0;
	var cl = [clipshadow(h0, glasses[gindex], object)]
	for (j=0;j<glasses.length;j++) {
		if (j !== gindex) {
			cl.push(clipshadow(h0, glasses[gindex], glasses[j]))
			cl.push(cliplight(h0, glasses[j], mu))
		}
	}
	return mergeandclip(glasses[gindex], cl);
}

function pointinrange(point, range) {
	for (var i=0;i<range.length;i++) {
		// xyarray = zip(range[i])
		if (pointInPolygon(point.x,point.y,range[i])) {
			return true;}
	}
	return false;
}

function PrepareData(walls, pC, pD, lw) {
	var pO = walls[0][0];
	var i = 0;
	var n = 0;
	var pA, pB;
	for (i=0;i<walls.length;i++) {
		// walls[i] = [pA, pB, Height,] + n
		pA = walls[i][0];
		pB = walls[i][1];
		pA.x = recx(pO.la,0,pA.la,0);
		pA.y = recy(0,pO.ln,pA.la,pA.ln);
		pB.x = recx(pO.la,0,pB.la,0);
		pB.y = recy(0,pO.ln,pB.la,pB.ln);
		n = 0;
		if (pA.la !== pB.la) {
			n = Math.acos(Math.asin(Math.sqrt(Math.pow((Math.sin(pB.la)-Math.sin(pA.la)),2)+Math.pow((Math.cos(pB.la)-Math.cos(pA.la)),2))/2)/Math.asin(Math.sqrt(Math.pow((Math.sin(pB.la)-Math.sin(pA.la)),2)+Math.pow((Math.cos(pB.la)-Math.cos(pA.la)*Math.cos(Math.abs(pA.ln-pB.ln))),2)+Math.pow((Math.cos(pA.la)*Math.sin(Math.abs(pA.ln-pB.ln))),2))/2));
		}
		if (pA.la<pB.la) {n=Math.PI-n}; if (pA.ln<pB.ln) {n=-n};
		walls[i].push(n);
	}
	pC.x = recx(pO.la,0,pC.la,0);
	pC.y = recy(0,pO.ln,pC.la,pC.ln);
	pD.x = recx(pO.la,0,pD.la,0);
	pD.y = recy(0,pO.ln,pD.la,pD.ln);
	var m = 0;
	if (pC.la !== pD.la) {
		m = Math.acos(Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)),2))/2)/Math.asin(Math.sqrt(Math.pow((Math.sin(pD.la)-Math.sin(pC.la)),2)+Math.pow((Math.cos(pD.la)-Math.cos(pC.la)*Math.cos(Math.abs(pC.ln-pD.ln))),2)+Math.pow((Math.cos(pC.la)*Math.sin(Math.abs(pC.ln-pD.ln))),2))/2));
	}
	if (pC.la<pD.la) {m=Math.PI-m}; if (pC.ln<pD.ln) {m=-m};
	// return pE, m
	return [{x: (pD.x*lw+pC.x*(100-lw))/100, y: (pD.y*lw+pC.y*(100-lw))/100}, m];
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

function Calc(daynum,time,pAvg,walls,pC,pD,pE,m,Hg,Hp,draw) {
var localtime = time+(pAvg.ln/pi180-120)/15;
var a = Math.PI*2*(daynum+localtime/24)/365;
var phi = 0.006918-0.399912*Math.cos(a)+0.070257*Math.sin(a)-0.006758*Math.cos(2*a)+0.000907*Math.sin(2*a)-0.002697*Math.cos(3*a)+0.001480*Math.sin(3*a);
var h0=Math.asin(Math.cos((localtime-12)*15*pi180)*Math.cos(pAvg.ln)*Math.cos(phi)+Math.sin(pAvg.ln)*Math.sin(phi));
if (h0 <= 0) {return -1;}
var cosA = (Math.sin(h0)*Math.sin(pAvg.la)-Math.sin(phi))/Math.cos(h0)/Math.cos(pAvg.la);
var alpha = Math.acos(cosA);
if (localtime < 12) {alpha = -alpha};
var n = walls[i][3];
var pA, pB, n, L, lightrange, worst = -1;
for (var i=0;i<walls.length;i++) {
	pA = walls[i][0];
	pB = walls[i][1];
	L = walls[i][2]/Math.tan(h0);
	n = walls[i][3];
	// 阳光照在背面
	if (Math.sin(alpha-n)>0) {if(draw !== -1) {return -2} else {return 0};}
	var mu = 2*n-alpha;
	lightrange = doclips(h0, i, glasses, [pC, pD, Hp], mu);
	if (draw !== -1) {lineGraph(draw, lightrange);}
	// 反射光照在背面
	if (Math.sin(mu-m)>0) {if (worst<0) {worst = 0}; continue;}
	if (pA.ln>pB.ln) {pB.y=-pB.y}; if (pA.la<pB.la) {pB.x=-pB.x};
	if (pA.ln>pC.ln) {pC.y=-pC.y}; if (pA.la<pC.la) {pC.x=-pC.x};
	if (pA.ln>pD.ln) {pD.y=-pD.y}; if (pA.la<pD.la) {pD.x=-pD.x};
	// 问：ela,eln,x5,y5 是什么？
	//答：ela，eln应该是窗的经纬度坐标，但是这里有问题，这里你用的公式好像是受照墙中点的公式。x4，y4应该是窗的直角坐标系坐标，但是不应该是由ela和eln算出来的，应该是用【附录p17-18 窗坐标公式】算出来的。x5，y5是照射到窗的那道反射光与产生这道反射光的玻璃幕墙建筑的交点。
	pF = {x: (pE.y-Math.tan(Math.PI-mu))/(pB.y/pB.x-Math.tan(Math.PI-mu)),
		  y: (pE.y-Math.tan(Math.PI-mu))*pB.y/((pB.y/pB.x-Math.tan(Math.PI-mu))*pB.x)};
	d = Math.sqrt(Math.pow(pE.x-pF.x,2)+Math.pow(pE.y-pF.y,2));
	hmax = Math.tan(h0)*(L-d);
	// 照不到高度
	if (hmax<Hp) {if (worst<0) {worst = 0}; continue};
	// 问：这是什么？和 pA2.la,pA2.ln,pB2.la,pB2.ln 有什么关系？
	//答：xa，ya和xb，yb分别是反射光投影终点的直角坐标系坐标【对应附录p17 反射光线投影终点坐标公式】；(draw !== -1) console.log([(pD.x*lw+pC.x*(100-lw))/100,(pD.y*lw+pC.y*(100-lw))/100,[0,xa,xb,pB.x],[0,ya,yb,pB.y]]);这段代码我看不懂，但是我知道公式是用来算窗的直角坐标系坐标的【附录p17-18 窗坐标公式】
	xa=L*Math.cos(Math.PI-mu);
	xb=pB.x+L*Math.cos(Math.PI-mu);
	ya=L*Math.sin(Math.PI-mu);
	yb=pB.y+L*Math.sin(Math.PI-mu);
	// if (draw !== -1) console.log([(pD.x*lw+pC.x*(100-lw))/100,(pD.y*lw+pC.y*(100-lw))/100,[0,xa,xb,pB.x],[0,ya,yb,pB.y]]);
	// 点不在范围
	// 问：此范围是之前算出来遮挡的总的范围？
	//答：范围是指受照射的范围吗？应该是算出一堵墙的反射光投影范围，然后算一遍窗是否在这个范围内，每堵墙都要算一遍。
	if (!pointinrange(pE, lightrange)) {if (worst<0) {worst = 0}; continue};

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
	if (worst<0) {worst = 0}; continue;
}
return worst;
}




function CalcOne(h0,alpha,pAvg,pA,pB,pA2,pB2,pC,pD,Hp,lw,draw,gindex,glasses) {
//daynum -> in one year
//time -> 24h, DEG
//ala,aln,bla,bln,cla,cln,dla,dln -> RAD!!!!!
// 新点坐标表示: {la, ln} or {la, ln, x, y}
// var pAvg = {la: (pA.la+pB.la)/2, ln: (pA.ln+pB.ln)/2};
// 确定范围
var lightrange = doclips(h0, gindex, glasses, )
if (draw !== -1) {
	pA2 = {la: pA.la-(L*Math.cos(Math.PI-mu)*2/dearth),
	       ln: pA.ln+(L*Math.sin(Math.PI-mu)*2/dearth)};
	pB2 = {la: pB.la-(L*Math.cos(Math.PI-mu)*2/dearth),
	       ln: pB.ln+(L*Math.sin(Math.PI-mu)*2/dearth)};
	lineGraph(draw,pA.la,pA.ln,pA2.la,pA2.ln,pB2.la,pB2.ln,pB.la,pB.ln);
}
// 反射光照在背面
if (Math.sin(mu-m)>0) {return 0;}
pB.x = recx(pA.la,0,pB.la,0);
pB.y = recy(0,pA.ln,pB.la,pB.ln);
pC.x = recx(pA.la,0,pC.la,0);
pC.y = recy(0,pA.ln,pC.la,pC.ln);
pD.x = recx(pA.la,0,pD.la,0);
pD.y = recy(0,pA.ln,pD.la,pD.ln);
if (pA.ln>pB.ln) {pB.y=-pB.y}; if (pA.la<pB.la) {pB.x=-pB.x};
if (pA.ln>pC.ln) {pC.y=-pC.y}; if (pA.la<pC.la) {pC.x=-pC.x};
if (pA.ln>pD.ln) {pD.y=-pD.y}; if (pA.la<pD.la) {pD.x=-pD.x};
// 问：ela,eln,x5,y5 是什么？
//答：ela，eln应该是窗的经纬度坐标，但是这里有问题，这里你用的公式好像是受照墙中点的公式。x4，y4应该是窗的直角坐标系坐标，但是不应该是由ela和eln算出来的，应该是用【附录p17-18 窗坐标公式】算出来的。x5，y5是照射到窗的那道反射光与产生这道反射光的玻璃幕墙建筑的交点。
pF = {x: (pE.y-Math.tan(Math.PI-mu))/(pB.y/pB.x-Math.tan(Math.PI-mu)),
      y: (pE.y-Math.tan(Math.PI-mu))*pB.y/((pB.y/pB.x-Math.tan(Math.PI-mu))*pB.x)};
d = Math.sqrt(Math.pow(pE.x-pF.x,2)+Math.pow(pE.y-pF.y,2));
hmax = Math.tan(h0)*(L-d);
// 照不到高度
if (hmax<Hp) {return 0};
// 问：这是什么？和 pA2.la,pA2.ln,pB2.la,pB2.ln 有什么关系？
//答：xa，ya和xb，yb分别是反射光投影终点的直角坐标系坐标【对应附录p17 反射光线投影终点坐标公式】；(draw !== -1) console.log([(pD.x*lw+pC.x*(100-lw))/100,(pD.y*lw+pC.y*(100-lw))/100,[0,xa,xb,pB.x],[0,ya,yb,pB.y]]);这段代码我看不懂，但是我知道公式是用来算窗的直角坐标系坐标的【附录p17-18 窗坐标公式】
xa=L*Math.cos(Math.PI-mu);
xb=pB.x+L*Math.cos(Math.PI-mu);
ya=L*Math.sin(Math.PI-mu);
yb=pB.y+L*Math.sin(Math.PI-mu);
// if (draw !== -1) console.log([(pD.x*lw+pC.x*(100-lw))/100,(pD.y*lw+pC.y*(100-lw))/100,[0,xa,xb,pB.x],[0,ya,yb,pB.y]]);
// 点不在范围
// 问：此范围是之前算出来遮挡的总的范围？
//答：范围是指受照射的范围吗？应该是算出一堵墙的反射光投影范围，然后算一遍窗是否在这个范围内，每堵墙都要算一遍。
if (!pointinrange(pE, ))

if (!pointInPolygon((pD.x*lw+pC.x*(100-lw))/100,(pD.y*lw+pC.y*(100-lw))/100,[0,xa,xb,pB.x],[0,ya,yb,pB.y])) {return 0};
gamma=-Math.PI/2+m-2*n+alpha;
lambda=Math.acos(Math.abs(Math.cos(gamma)*Math.cos(h0)));
B=rho*137000*Math.exp(-0.223/Math.sin(h0))/Math.PI;

if (lambda<Math.PI/12) {
	if (B>=2000) {return 5;};
} else if (lambda>=Math.PI/12 && lambda<=Math.PI/6) {
	if (B>=2000 && B<4000) {
		return 2;
	} else if (B>=4000 && B<6000) {
		return 3;
	} else if (B>=6000) {
		return 4;
	};
} else if (lambda>Math.PI/6) {
	return 1;
};
return 0;
}

