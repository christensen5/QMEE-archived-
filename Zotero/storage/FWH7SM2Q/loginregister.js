$(document).ready(function () {
    
	// When radio changes
    $('input[name="accountTypeLogin[type]"]').click(function () {  
        var val = $('input[name="accountTypeLogin[type]"]:checked').val();
        if(val == '2') {
        	// Change to Corp Form    	
        	makeCorpLogin();
        } else {
        	// RingGo form
        	makeRingGoLogin();
        }
    });
    
    
    $('input[name="accountType[type]"]').click(function () {  
        var val = $('input[name="accountType[type]"]:checked').val();
        if(val == '2') {
        	$('.permitRadio').next().addClass('permitRegister');
                $('.parkingRadio').next().removeClass('parkingRegister');
        } else {
        	$('.permitRadio').next().removeClass('permitRegister');
                $('.parkingRadio').next().addClass('parkingRegister');
        }
    });
    
    // Check current form on pageload
    var val = $('input[name="accountTypeLogin[type]"]:checked').val();
    if(val == '2') {
    	//A fix for IE
	$('.corporateRadio').next().addClass('corporateLogin');
	$('.personalRadio').next().removeClass('personalLogin');
        
        $('.permitRadio').next().addClass('permitRegister');
        $('.parkingRadio').next().removeClass('parkingRegister');
    } else {
    	//A fix for IE
	$('.corporateRadio').next().removeClass('corporateLogin');
	$('.personalRadio').next().addClass('personalLogin');
        
        $('.permitRadio').next().removeClass('permitRegister');
        $('.parkingRadio').next().addClass('parkingRegister');
    }
    
    $('form#Login').attr('onsubmit','checkSubmitType()');
    $('#currentMembers input[type="submit"]').attr('onclick','submittedButton=1;');
    $('#newMembers input[type="submit"]').attr('onclick','submittedButton=2;');
    
    //Handler for WCC new Log In page manipulations.
    if(url.toLowerCase().indexOf("westminster") >= 0) {
        $('body.wcc-new').addClass('wcc-login-page');
        $('.wcc-new.wcc-login-page input[name="labyrinth_Login_next"]').parent().attr("style", "text-align: left; vertical-align: text-top;");
        $('.wcc-new.wcc-login-page input#field-pin').attr("placeholder", "Password");
        $('.wcc-new.wcc-login-page input#field-cli').attr("placeholder", "e.g. jbloggs@gmail.com");
        $('.wcc-new.wcc-login-page div.element.div-static1 > div').html('<span class="wcc-forgot-text"><a href="/resetpassword">Forgotten password or PIN?</a>');
        $('.wcc-new.wcc-login-page div.element.div-static2 > div').html('First time? <a href="/register">Register here</a>');
        $('.wcc-new.wcc-login-page .tablelessform li span.error').find('br').remove();
        $('.wcc-new.wcc-login-page li.cli-wrapper label#label-cli').html("Email or phone:");
        $('.wcc-new.wcc-login-page li.pin-wrapper label#label-pin').html("Password, PIN:");
        if ($('.wcc-new.wcc-login-page .warning-notification > p') > -1) {
                var headerErrorText = $('.wcc-new.wcc-login-page .warning-notification > p')[0].html();
                $('.wcc-new.wcc-login-page .warning-notification').attr("style", "display: none;");
                $('.wcc-new.wcc-login-page .accountTypeLogin-wrapper').after('<li class="static" style="text-align: center; font-weight: bold;"><div class="element"><div class="elemErrorWrap" style="color: #d4463f;">' + headerErrorText + '</div></div></li>');
        }

        function validateEmail(obj) {
            removeValidityIcons($(obj));
            if(/^(([\w\.]*)\@(\S*\.)?\w+)$|^(\d*)$/.test($(obj).val()) && $(obj).val() != '') {
                $(obj).after('<div class="position-validity"><img class="valid" src="/images/westminster/ico-tick-green.png" /></div>');
            } else {
                $(obj).after('<div class="position-validity"><img class="invalid" src="/images/westminster/ico-cross-red.svg" /></div>');
            }
        }

        function validatePassword(obj) {
            removeValidityIcons($(obj));
            if($(obj).val() != '') {
                $(obj).after('<div class="position-validity"><img class="valid" src="/images/westminster/ico-tick-green.png" /></div>');
            } else {
                $(obj).after('<div class="position-validity"><img class="invalid" src="/images/westminster/ico-cross-red.svg" /></div>');
            }
        };

        function removeErrors(obj) {
            $(obj).closest('.error').removeClass('error');
        }

        function removeValidityIcons(obj) {
            $(obj).next('.position-validity').remove();
        };

        $('.wcc-new.wcc-login-page li.cli-wrapper input#field-cli').blur(function() {
            validateEmail($(this));
        })

        $('.wcc-new.wcc-login-page li.cli-wrapper input#field-cli').focus(function(){
            removeErrors($(this));
        });

        $('.wcc-new.wcc-login-page li.pin-wrapper input#field-pin').blur(function(){
            validatePassword($(this));
        })

        $('.wcc-new.wcc-login-page li.pin-wrapper input#field-pin').focus(function(){
            removeErrors($(this));
            validatePassword($(this));
        });

        // setTimeout(function(){
        //     $('li.cli-wrapper input#field-cli').keyup();
        // }, 300);
    }
});

// Store default lab action & get URL for Corp Login
var currentAction = $('form#Login').attr('action');
var url = document.URL;
var CorpURL;
if(url.toLowerCase().indexOf("westminster") >= 0) {
	CorpURL = window.location.protocol+'//westminster-corporate.'+document.location.hostname+'/corplogin';
} else {
	CorpURL = window.location.protocol+'//corporate.'+document.location.hostname+'/corplogin';
}
var submittedButton; // 1 is RingGo, 2 is Corp

// Check Submission Type when form clicked and set action as appropriate
function checkSubmitType() {
	if(submittedButton == 1) {
		var val = $('input[name="accountTypeLogin[type]"]:checked').val();
		
		if(val == '2') {
			// Change to Corp Form
	    	$('form#Login').attr('action',CorpURL);
		} else {
			// RingGo form
			$('form#Login').attr('action',currentAction);
		}
	}
}

function makeCorpLogin() {
	$('#field-cli').attr('name', 'Email');
	$('#field-pin').attr('name', 'Password');
	
	$('#label-cli').text('Email address');
	$('#label-pin').text('Password');
	
	$('.personalRadio').next().removeClass('personalLogin');
	$('.corporateRadio').next().addClass('corporateLogin');
}

function makeRingGoLogin() {
	$('#field-cli').attr('name', 'cli');
	$('#field-pin').attr('name', 'pin');
	
	$('#label-cli').text('Mobile number or E-mail');
	$('#label-pin').text('Account PIN/Password');
	
	//A fix for IE
	$('.corporateRadio').next().removeClass('corporateLogin');
	$('.personalRadio').next().addClass('personalLogin');
}