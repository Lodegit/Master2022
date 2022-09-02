// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, true, false, false, false, false, false, false, false, true, true, false, false, true ];
var arrayMetadata    = [ [ "1", "GSM1244762_BH2010610-8-8PCyuan.CEL", "GSM1244762", "adjacent_tissue", "P09", "10/15/11 16:33:01" ], [ "2", "GSM1244764_BH2010610-8-8TCyuan.CEL", "GSM1244764", "carcinoma", "P09", "10/15/11 16:44:30" ], [ "3", "GSM1244766_BH2010610-8-10TCyuan.CEL", "GSM1244766", "carcinoma", "P10", "10/15/11 11:39:42" ], [ "4", "GSM1244768_BH2010610-8-10PCyuan.CEL", "GSM1244768", "adjacent_tissue", "P10", "10/19/11 11:32:20" ], [ "5", "GSM1244769_BH2010610-8-16TCyuan.CEL", "GSM1244769", "carcinoma", "P11", "10/15/11 13:25:32" ], [ "6", "GSM1244770_BH2010610-8-16PCyuan.CEL", "GSM1244770", "adjacent_tissue", "P11", "10/15/11 11:28:35" ], [ "7", "GSM1244775_BH2010610-8-18PCyuan.CEL", "GSM1244775", "adjacent_tissue", "P12", "10/19/11 12:59:32" ], [ "8", "GSM1244776_BH2010610-8-18TCyuan.CEL", "GSM1244776", "carcinoma", "P12", "10/19/11 13:21:53" ], [ "9", "GSM1244777_BH2010610-8-20TCyuan.CEL", "GSM1244777", "carcinoma", "P13", "10/19/11 14:56:47" ], [ "10", "GSM1244778_BH2010610-8-20PCyuan.CEL", "GSM1244778", "adjacent_tissue", "P13", "10/19/11 14:43:46" ], [ "11", "GSM1244782_BH2010610-8-46TC.CEL", "GSM1244782", "carcinoma", "P14", "10/14/11 10:30:04" ], [ "12", "GSM1244783_BH2010610-8-46PC.CEL", "GSM1244783", "adjacent_tissue", "P14", "10/14/11 10:52:27" ], [ "13", "GSM1244785_BH2010610-8-58TC.CEL", "GSM1244785", "carcinoma", "P15", "10/14/11 12:35:44" ], [ "14", "GSM1244787_BH2010610-8-58PC.CEL", "GSM1244787", "adjacent_tissue", "P15", "10/19/11 16:19:56" ], [ "15", "GSM1244790_BH2010610-8-62TC.CEL", "GSM1244790", "carcinoma", "P16", "10/14/11 16:32:19" ], [ "16", "GSM1244792_BH2010610-8-62PC.CEL", "GSM1244792", "adjacent_tissue", "P16", "10/14/11 13:50:38" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
