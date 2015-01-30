/*
{
	var nav1off = new Image();
	nav1off.src = "en/images/nav_home_off.gif";
	var nav1on = new Image();
	nav1on.src = "en/images/nav_home_on.gif";
	var nav2off = new Image();
	nav2off.src = "en/images/nav_guide_off.gif";
	var nav2on = new Image();
	nav2on.src = "en/images/nav_guide_on.gif";
	var nav3off = new Image();
	nav3off.src = "en/images/nav_loffice_off.gif";
	var nav3on = new Image();
	nav3on.src = "en/images/nav_loffice_on.gif";
	var nav4off = new Image();
	nav4off.src = "en/images/nav_soffice_off.gif";
	var nav4on = new Image();
	nav4on.src = "en/images/nav_soffice_on.gif";
	var nav5off = new Image();
	nav5off.src = "en/images/nav_research_off.gif";
	var nav5on = new Image();
	nav5on.src = "en/images/nav_research_on.gif";
	var nav6off = new Image();
	nav6off.src = "en/images/nav_sdist_off.gif";
	var nav6on = new Image();
	nav6on.src = "en/images/nav_sdist_on.gif";
	var nav7off = new Image();
	nav7off.src = "en/images/nav_ldist_off.gif";
	var nav7on = new Image();
	nav7on.src = "en/images/nav_ldist_on.gif";
	var nav8off = new Image();
	nav8off.src = "en/images/nav_manu_off.gif";
	var nav8on = new Image();
	nav8on.src = "en/images/nav_manu_on.gif";
	var nav9off = new Image();
	nav9off.src = "en/images/nav_glossary_off.gif";
	var nav9on = new Image();
	nav9on.src = "en/images/nav_glossary_on.gif";
	var nav10off = new Image();
	nav10off.src = "en/images/nav_index_off.gif";
	var nav10on = new Image();
	nav10on.src = "en/images/nav_index_on.gif";
	var nav11off = new Image();
	nav11off.src = "en/images/nav_tools_off.gif";
	var nav11on = new Image();
	nav11on.src = "en/images/nav_tools_on.gif";

	var hombody0off = new Image();
	hombody0off.src = "en/images/hom_largeroffices_off.gif";
	var hombody0on = new Image();
	hombody0on.src = "en/images/hom_largeroffices_on.gif";
	var hombody1off = new Image();
	hombody1off.src = "en/images/hom_smalloffices_off.gif";
	var hombody1on = new Image();
	hombody1on.src = "en/images/hom_smalloffices_on.gif";
	var hombody2off = new Image();
	hombody2off.src = "en/images/hom_research_off.gif";
	var hombody2on = new Image();
	hombody2on.src = "en/images/hom_research_on.gif";
	var hombody3off = new Image();
	hombody3off.src = "en/images/hom_sdistribution_off.gif";
	var hombody3on = new Image();
	hombody3on.src = "en/images/hom_sdistribution_on.gif";
	var hombody4off = new Image();
	hombody4off.src = "en/images/hom_ldistribution_off.gif";
	var hombody4on = new Image();
	hombody4on.src = "en/images/hom_ldistribution_on.gif";
	var hombody5off = new Image();
	hombody5off.src = "en/images/hom_manufacturing_off.gif";
	var hombody5on = new Image();
	hombody5on.src = "en/images/hom_manufacturing_on.gif";

}
*/

var isNav6, isNav4, isIE;

if (document.getElementById) {
	isNav6 = true;
}

else if (document.layers) {
	isNav4 = true;
}

else {
	isIE = true;
}

/* check browser type
 * supported  - Netscape 4.5 or above
 *            - Internet Explorer 4.0 or above
 */
function checkBrowser ()
{
  // declare variables
  var strBrowser  = navigator.appName;
  var fltVer = parseFloat(navigator.appVersion);
  var agt=navigator.userAgent.toLowerCase();
  var is_ie     = ((agt.indexOf("msie") != -1) && (agt.indexOf("opera") == -1));
  var is_ie4up  = (is_ie && (fltVer >= 4));
  
  var is_nav  = ((agt.indexOf('mozilla')!=-1) && (agt.indexOf('spoofer')==-1)
                && (agt.indexOf('compatible') == -1) && (agt.indexOf('opera')==-1)
                && (agt.indexOf('webtv')==-1) && (agt.indexOf('hotjava')==-1));
  var is_nav45up = (is_nav && (fltVer >= 4.5));

  // determine content based on the browswer being used
  //alert("Your browser is: " + strBrowser + " " + fltVer);
  if ((strBrowser == "Netscape" && !is_nav45up) 
    || (strBrowser == "Microsoft Internet Explorer" && !is_ie4up))
  {
    spawnPDFWindow('/site/browser_check.affx');
  }
}

function changeImage(imgDocID,imgObjName) {
	if (document.images[imgDocID].complete) {
		document.images[imgDocID].src = eval(imgObjName + ".src");
	}
}

function fixIt() {
	if (isNav4) {
		origWidth = innerWidth;
		origHeight = innerHeight;
		if (innerWidth != origWidth || innerHeight != origHeight) {
			location.reload();
		}
	}
}

function arrayfinderPopup() {
  window.open('array_finder.html','ArrayFinder','width=500,height=400,screenX=100,screenY=100,top=100,left=100,resizable=yes,toolbar=no,scrollbars=yes');
}

function enlarge(path) {
  window.open ("/community/wayahead/enlarge.jsp?image="+path, "Affymetrix_"+which, "height=500,width=500");
  which++;
}

function enlargeImg() {
	var query = window.location.search.substring(1);
	var values = query.split("=");
	var path = values[1]; //returns image path
	var stringA = path.substring(0,(path.length-4));
	var ext = path.substring((path.length-4),path.length);
	var newImg = stringA + "_enlarge" + ext;
	return newImg;
}

function spawnWindow(url, windowname, attributes, mainwindowname) {
  remote = window.open(null,windowname, attributes);
  remote.location = url;
  remote.focus();
  if (remote.opener == null) {
    remote.opener = window;
  }
  if (remote.opener.name == null){
  remote.opener.name = mainwindowname;
  }
}

function spawnHelpWindow(url) {
  help = spawnWindow(
    url,
    "affx_help",
    "width=500,height=550,screenX=100,screenY=100,top=100,left=100,resizable=yes,toolbar=no,scrollbars=yes",
    "affymain");
}

function spawnPDFWindow(url) {
  help = spawnWindow(
    url,
    "affx_pdf",
    "width=750,height=550,screenX=150,screenY=150,top=150,left=150,resizable=yes,toolbar=yes,location=yes,scrollbars=yes",
    "");
}
