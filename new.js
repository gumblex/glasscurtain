function recx(la0,_ln0,lax,_lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.sin(lax)-Math.sin(la0),2)+Math.pow(Math.cos(lax)-Math.cos(la0),2))/2)};
function recy(_la0,ln0,lax,lnx) {return dearth*Math.asin(Math.sqrt(Math.pow(Math.cos(lax)-Math.cos(lax)*Math.cos(Math.abs(ln0-lnx)),2)+Math.pow(Math.cos(lax)*Math.sin(Math.abs(ln0-lnx)),2))/2)};

function checkclockwise(points){
	var sum = 0;
	for (var i=1;i<points.length;i++){
		sum += (points[i][0]-points[i-1][0])*(points[i][1]+points[i-1][1]);}
	sum += (points[i-1][0]-points[0][0])*(points[i-1][1]+points[0][1])
	if (sum < 0){return points.reverse()} else {return points};
}


// y = Math.tan(-alpha) * (x - x4) + y4;
// y = Math.tan(-alpha) * (x - x5) + y5;
// y = (y7 - y6) / (x7 - x6) * (x - x6) + y6;
function clipshadow(){
	x8 = (y6 - y4 - (x6 * (y7 - y6) / (x7 - x6)) + Math.tan(-alpha) * x4) / (Math.tan(-alpha) - (y7 - y6) / (x7 - x6));
	y8 = Math.tan(-alpha) * (y6 - y4 - (x6 * (y7 - y6) / (x7 - x6)) + Math.tan(-alpha) * x4) / (Math.tan(-alpha) - (y7 - y6) / (x7 - x6)) + y4;
	x9 = (y6 - y5 - (x6 * (y7 - y6) / (x7 - x6)) + Math.tan(-alpha) * x5) / (Math.tan(-alpha) - (y7 - y6) / (x7 - x6));
	y9 = (Math.tan(-alpha) * (y6 - y5 - (x6 * (y7 - y6) / (x7 - x6)) + Math.tan(-alpha) * x5)) / (Math.tan(-alpha) - (y7 - y6) / (x7 - x6)) + y5;
	if ((x8 < x6 || x8 > x7) && (x9 < x6 || x9 > x7)) {return 0;}

	d1 = Math.sqrt(Math.pow(x4 - x8, 2) + Math.pow(y4 - y8, 2));
	d2 = Math.sqrt(Math.pow(x5 - x9, 2) + Math.pow(y5 - y9, 2));
	L4 = H / Math.tan(ho);
	Hmax1 = Math.tan(ho) * (L4 - d1);
	Hmax2 = Math.tan(ho) * (L4 - d2);
	if (Hmax1 < 0 && Hmax2 < 0) {return 0;}
	L1 = Hmax1 / Math.tan(ho);
	L2 = Hmax2 / Math.tan(ho);
	la5 = laE - (d1 * Math.cos(Math.PI - alpha)) / r;
	ln5 = lnE + (d1 * Math.sin(Math.PI - alpha)) / r;
	la6 = laF - (d2 * Math.cos(Math.PI - alpha)) / r;
	ln6 = lnF + (d2 * Math.sin(Math.PI - alpha)) / r;
	la7 = la5 - (L1 * Math.cos(Math.PI - mu)) / r;
	ln7 = ln5 + (L1 * Math.sin(Math.PI - mu)) / r;
	la8 = la6 - (L2 * Math.cos(Math.PI - mu)) / r;
	ln8 = ln6 + (L2 * Math.sin(Math.PI - mu)) / r;
	if (Hmax1 > 0 && Hmax2 > 0) {
		range1 = clip_polygon_difference([[ala,aln],[ala2,aln2],[bla2,bln2],[bla,bln]], [[la5,ln5],[la6,ln6],[la8,ln8],[la7,ln7]]);
	} else if (Hmax1 < 0) {
		range1 = clip_polygon_difference([[ala,aln],[ala2,aln2],[bla2,bln2],[bla,bln]], [[la6,ln6],[la8,ln8],[la7,ln7]]);
	} else {
		range1 = clip_polygon_difference([[ala,aln],[ala2,aln2],[bla2,bln2],[bla,bln]], [[la5,ln5],[la8,ln8],[la7,ln7]]);
	}
	L3 = H / Math.tan(ho)
}
function cliplight(){
	la11 = la9 - (L1 * Math.cos(Math.PI - mu)) / r
	ln11 = ln9 + (L1 * Math.sin(Math.PI - mu)) / r
	la12 = la10 - (L2 * Math.cos(Math.PI - mu)) / r
	ln12 = ln10 + (L2 * Math.sin(Math.PI - mu)) / r
	var range2 = [];
	for (var k=0;k<range1.length;k++) {
		range2.push(clip_polygon_difference(range1[k], [[la9,ln9],[la10,ln10],[la12,ln12],[la11,ln11]]));
	}
}

function doclips(glasses, object){
	
}



function clip (subjectPolygon, clipPolygon) {

	var cp1, cp2, s, e;
	var inside = function (p) {
		return (cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0]);
	};
	var intersection = function () {
		var dc = [cp1[0] - cp2[0], cp1[1] - cp2[1]],
			dp = [s[0] - e[0], s[1] - e[1]],
			n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0],
			n2 = s[0] * e[1] - s[1] * e[0], 
			n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0]);
		return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3];
	};
	var outputList = subjectPolygon;
	cp1 = clipPolygon[clipPolygon.length-1];
	for (j in clipPolygon) {
		var cp2 = clipPolygon[j];
		var inputList = outputList;
		outputList = [];
		s = inputList[inputList.length - 1]; //last on the input list
		for (i in inputList) {
			var e = inputList[i];
			if (inside(e)) {
				if (!inside(s)) {
					outputList.push(intersection());
				}
				outputList.push(e);
			}
			else if (inside(s)) {
				outputList.push(intersection());
			}
			s = e;
		}
		cp1 = cp2;
	}
	return outputList
}

function drawPolygon(context, polygon, strokeStyle, fillStyle) {
	context.strokeStyle = strokeStyle;
	context.fillStyle = fillStyle;
	context.beginPath();
	context.moveTo(polygon[0][0],polygon[0][1]); //first vertex
	for (var i = 1; i < polygon.length ; i++)
		context.lineTo(polygon[i][0],polygon[i][1]);
	context.lineTo(polygon[0][0],polygon[0][1]); //back to start
	context.fill();
	context.stroke();
	context.closePath();
}

window.onload = function () {
	var context = document.getElementById('canvas').getContext('2d');
	var subjectPolygon = [[50, 150], [200, 50], [350, 150], [350, 300], [250, 300], [200, 250], [150, 350], [100, 250], [100, 200]],
		clipPolygon = [[100, 100], [300, 100], [300, 300], [100, 300]];
	var clippedPolygon = clip(subjectPolygon, clipPolygon);
	drawPolygon(context, clipPolygon, '#888','#88f');
	drawPolygon(context, subjectPolygon, '#888','#8f8');
	drawPolygon(context, clippedPolygon, '#000','#0ff');
}
