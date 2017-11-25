_satellite.pushBlockingScript(function(event, target, $variables){
  // milliseconds to midnight function
// useful for using _satellite.setCookie() to expire at 11:59 PM tonight
// _satellite.setCookie('cname','cvalue', _sdi.msToMidnight / (24*60*60*1000)

_sdi.msToMidnight = function() {
	var now = new Date();
	var hours = (24-now.getHours())*60*60;
	var minutes = (60-now.getMinutes())*60;
	var seconds = 60-now.getSeconds();
	var milliseconds = ( hours + minutes + seconds ) * 1000; 
  
  return milliseconds;
}
});
