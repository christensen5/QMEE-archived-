_satellite.pushBlockingScript(function(event, target, $variables){
  // return a three character representation of today's date

_sdi.dateCode3 = function( d ) {
  var today = d ? new Date(d) : new Date();
	var monthCodes = '1234567890AB';
	var dayCodes = '1234567890ABCDEFGHIJKLMNOPQRSTU';
	return( 
         today.getFullYear().toString().charAt(3) 
         + monthCodes.charAt(today.getMonth()) 
         + dayCodes.charAt(today.getDate()-1) 
        );
}
});
