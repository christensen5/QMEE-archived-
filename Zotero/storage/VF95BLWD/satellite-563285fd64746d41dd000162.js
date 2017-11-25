//console.log("oas.js: oas_tag.query = "+oas_tag.query);

function hawkeye_query(target) {
    var result = '';
    var content = target.content;
    var content_length = content.length;
    for (var i = 0; i < content_length; i++) {
       var item = content[i];
       if (item != '') {
          if (result != '') {
              result += '&';
          }
          result += "content=" + item;
       }
    }
    return result;
}

function hawkeye_updateOASTagQuery(oasquery, target) {
    var queryaddon = '';
    var result = oasquery;
    if (typeof target === 'object') {
         queryaddon = hawkeye_query(target);
        if ((typeof oasquery !== 'undefined') && (oasquery != '')) {
            if (queryaddon != '') {
                result = oasquery + "&" + queryaddon;
            }
        } else {
            if (queryaddon != '') {
                result = queryaddon;
            }
        }
    } 
    return result;
}

if (typeof hawkeye_target !== 'undefined') {
	oas_tag.query = hawkeye_updateOASTagQuery(oas_tag.query, hawkeye_target);
}

//console.log("oas.js: oas_tag.query = "+oas_tag.query);

oas_tag.url       = 'oascentral.sciencemag.org'; // Define OAS URL

var site_page = _satellite.getVar('PageName_CustomScript').split("/")[0].replace('.:','.sciencemag.org');
var site_path = document.location.pathname;

if(location.host.indexOf('www.') === 0){
    site_page = location.host.replace('www.','');
}

var devsites = ["sciencedev.aaas.org",
"rjensen.sciencemag.org",
"scienceqa.aws.aaas.org",
"dtmdemoads.dev",
"sciencecore.dev",
"dmanalytics.aaas.org",
"sciblogsdev.aaas.org",
"dmdevweb01.aws.aaas.org",
"hw-f5-immunology.highwire.org",
"jnl-aaas-sci.drupal-stage-jnl-web02.highwire.org",
"sci-stage-aaascodev.stage.highwire.org"];

for (var i = 0; i < devsites.length; i++) {
    if (document.location.hostname === devsites[i]){
      site_page = "sciencedev.aaas.org";
      site_path = "";
    }
}

//add front if landing page
if(document.location.pathname === "/") {
	site_page += "/front";
	site_path = "";
}

//login pages
if(document.location.pathname === "/user/login") {
	site_page += "/front";
	site_path = "";
}

//GHP
if(document.location.hostname === "www.sciencemag.org" && document.location.pathname === "/") {
	site_page = "sciencemag.org";
	site_path = "/ghp"
}

//GHP (Dev)
if(document.location.hostname === "sciencedev.aaas.org" && document.location.pathname === "/") {
	site_page = "sciencedev.aaas.org";
	site_path = "/ghp"
}

if(site_page === "sciencemag.org" && document.location.pathname === "/news") {
 	site_page = "sciencemag.org/newshomepage";
 	site_path = "";
}

if(site_page === "sciencemag.org" && document.location.pathname === "/careers") {
 	site_page = "sciencemag.org/careershomepage";
 	site_path = "";
}

if(site_page === "blogs.sciencemag.org" && document.location.pathname === "/pipeline/") {
 	site_page = "blogs.sciencemag.org/pipelinehome";
 	site_path = "";
}

/*if(site_path.match(/\/features/)) {
 	site_path = "/librarian"; //force house ads
}*/

oas_tag.site_page = site_page + site_path;

if (oas_tag.site_page.match(/^\blocalhost{1}/)) {
  var local_page = oas_tag.site_page.replace(/^\blocalhost{1}/, "sciencedev.aaas.org");
  oas_tag.site_page = local_page;
  //oas_tag.reloadAds();
};

oas_tag.allowSizeOverride = true; // Needed for responsive sizes
oas_tag.lazyload = false; // Async loading
oas_tag.anchorElementPrefix = "aaas-oas"; // This customization has been requested from HW for Jcore
oas_tag.analytics = true;

(function() {
  oas_tag.version ='1'; oas_tag.loadAd = oas_tag.loadAd || function(){}; var oas = document.createElement('script'),
  protocol = 'https:' == document.location.protocol?'https://':'http://', node = document.getElementsByTagName('script')[0];
  oas.type = 'text/javascript'; oas.async = true;
  oas.src = protocol + oas_tag.url + '/om/' + oas_tag.version + '.js'; node.parentNode.insertBefore(oas, node);
})();

oas_tag.sizes = function() {
	var adCatch = document.querySelectorAll('[id^="aaas-oas_"]');
	var adActive = [];
	var i;
	for(i = 0; i < adCatch.length; i++){
	  adActive.push(adCatch[i].id);
	}
	var clip = adActive.toString();
	var resorter = clip.replace(/aaas-oas_/g,"");
	var regather = resorter.split(',');

	var j;
	var adWidth = 0;
	var adHeight = 0;

	for(j = 0; j < regather.length; j++){
	  if(regather[j] === "Top"){
		adWidth = 728;
		adHeight = 90;
	  }
	  else if(regather[j] === "Top2") {
		adWidth = 970;
		adHeight = 1;
	  }
	  else if(regather[j] === "Right2") {
		adWidth = 300;
		adHeight = 600;
	  }
	  else if(regather[j] === "x30") {
		adWidth = 1280;
		adHeight = 60;
	  }
	  else if(regather[j] === "x51") {
		adWidth = 300;
		adHeight = 250;
	  }
	  else if(regather[j] === "x60") {
		adWidth = 300;
		adHeight = 250;
	  }
	  else if(regather[j] === "x31") {
		adWidth = 100;
		adHeight = 100;
	  }
	  else{
	  adWidth = 300;
		adHeight = 125;
	  }
	  oas_tag.definePos(regather[j],	[adWidth, adHeight]);
	}
};

var sheet = (function() {
  var style = document.createElement("style");
  style.appendChild(document.createTextNode(""));
  document.head.appendChild(style);
  return style.sheet;
})();

sheet.insertRule("[id^=aaas-oas]:not(#aaas-oas_x30) { position: relative; }", 0);

sheet.insertRule(".ad__label {position: relative; font-size: .75rem; font-family: Benton Sans, Helvetica, Arial, sans-serif; margin: 0; text-align: center; line-height: 1.1; color: gray; }", 0);

sheet.insertRule("#aaas-oas_Top .ad__label {color: white;}", 0);

sheet.insertRule("#aaas-oas_Top2 .ad__label {padding: .2rem; text-align: left; top: 0; left: 0; position: absolute; color: gray; width: 25px;}", 0);

sheet.insertRule("#aaas-oas_Top2 { background-color: #ffffff; width: 970px; text-align: center; }", 0);

//sheet.insertRule("#aaas-oas_Top2 { background-color: #e6e6e6; width: 970px; text-align: center; position: static; }", 0);

//sheet.insertRule("#aaas-oas_Top2 .billboard--closead {font-size: .75rem; font-family: Helvetica, Arial, sans-serif; color: white; width: 100px; float: right; position: relative; top: 0; z-index: 2; margin-bottom: -20px; background-color: rgba(0, 0, 0, 0.5); cursor: pointer;}", 0);

sheet.insertRule("#aaas-oas_x30 .ad__label {display: inline-block;  top: 10; left: 0;  text-align: left; }", 0);

sheet.insertRule(".oas-position-x30 .container--footer {padding-bottom: 70px;}", 0);

sheet.insertRule("#aaas-oas_x51, #aaas-oas_Right2, #aaas-oas_x60 {text-align: center;}", 0);

sheet.insertRule("#aaas-oas_x30 {max-width: 1280px; line-height: 0; bottom: 0;  position: fixed; right: 0;  left: 0;  margin: 0px auto;  z-index: 10000;} ", 0);

sheet.insertRule("#aaas-oas_x30 a {line-height: 0;  display: block;} ", 0);

sheet.insertRule("#aaas-oas_x30 .ad__zapper {right: 10px; position: absolute; z-index: 1; cursor: pointer; } ", 0);

sheet.insertRule("#aaas-oas_x31 .ad__label {display: none; visibility: hidden;} ", 0);

sheet.insertRule("#aaas-oas_x31 {z-index: -100; width: 100px; height: 100px;} ", 0);

sheet.insertRule("@media (min-width: 961px) {#aaas-oas_Top .ad__label {display: inline-block; transform: rotate(-90deg); text-align: left; transform-origin: left bottom; bottom: 4px; left: 0; position: absolute; color: white; }}", 0);

sheet.insertRule("@media (min-width: 1279px) {#aaas-oas_x30 {width: 1280px !important; max-width: 1280px;} #aaas-oas_x30 .ad__zapper{right: -10px;}}", 0);

oas_tag.callbackHandler = function() {

  Element.prototype.remove = function() {
    this.parentElement.removeChild(this);
  }

  NodeList.prototype.remove = HTMLCollection.prototype.remove = function() {
    for(var i = this.length - 1; i >= 0; i--) {
      if(this[i] && this[i].parentElement) {
        this[i].parentElement.removeChild(this[i]);
      }
    }
  }

  oas_tag.addHandler('callbackId', function(data) {

    var pageAds = document.querySelectorAll("[id^=" + oas_tag.anchorElementPrefix + "]"),
        adPositions = {};

    for (var i=0; i < pageAds.length; i++) {
      posKey = pageAds[i].id;
      posVal = pageAds[i].id.split('_')[1];
      adPositions[posKey] = posVal;
    }

    for (var key in adPositions) {
      var posData = adPositions[key];

      document.body.classList.add("oas-position-" + adPositions[key]);

      if (data[posData] && data[posData].creativeId != "empty.gif") {
        if(data[posData].position === "Top2") {        
        	document.getElementById("aaas-oas_Top2").style.marginTop="15px";
        	document.getElementById("aaas-oas_Top2").style.marginBottom="15px";
        	document.getElementById("aaas-oas_Top2").style.maxHeight="426px";            
        }
        document.getElementById(key).insertAdjacentHTML("afterbegin", '<p class="ad__label">Advertisement</p>');
      }
      if (data[posData].position === "x30") {
        document.getElementById(key).insertAdjacentHTML("afterbegin", '<img onclick="this.parentNode.remove()" width="24" height="24" class="ad__zapper" src="data:image/svg+xml;utf8,<svg viewBox=\'0 0 22 22\' xmlns=\'http://www.w3.org/2000/svg\'><g transform=\'translate(1 1)\' stroke=\'#fff\' stroke-width=\'2\' fill=\'none\' fill-rule=\'evenodd\'><circle fill-opacity=\'.8\' fill=\'#000\' cx=\'10\' cy=\'10\' r=\'10\'/><g stroke-linecap=\'square\'><path d=\'M5.5 5.5l9 9M14.5 5.5l-9 9\'/></g></g></svg>">');
      }
    }
  })
};
