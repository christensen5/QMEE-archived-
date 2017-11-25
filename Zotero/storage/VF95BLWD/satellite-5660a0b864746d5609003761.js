_satellite.pushAsyncScript(function(event, target, $variables){
  insertChunk("//sciencestatic.aws.aaas.org.s3.amazonaws.com/dtm/form-subscribe-capture.html", "#aaas-util-foot-1");

var sub = "science";

var subdomain = document.location.hostname.match(/([a-zA-Z+]+)\./)[1];

if(subdomain === "advances" || subdomain === "stke" || subdomain === "stm" || subdomain === "immunology" || subdomain === "robotics") {
	sub = subdomain;
}

if(document.location.pathname === "/journals/robotics") {
	sub = "robotics";
}

insertChunk("//sciencestatic.aws.aaas.org.s3.amazonaws.com/dtm/footer-newsletter-form/form-" + sub + ".html", "#aaas-util-foot-2");
});
