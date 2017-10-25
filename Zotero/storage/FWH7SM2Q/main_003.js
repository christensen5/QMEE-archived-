//  For old browsers
(function () {
    var method;
    var noop = function () {};
    var methods = [
        'assert', 'clear', 'count', 'debug', 'dir', 'dirxml', 'error', 'exception', 'group', 'groupCollapsed', 'groupEnd', 'info', 'log', 'markTimeline', 'profile', 'profileEnd', 'table', 'time', 'timeEnd', 'timeline', 'timelineEnd', 'timeStamp', 'trace', 'warn'
    ];
    var length = methods.length;
    var console = (window.console = window.console || {});
    while (length--) {
        method = methods[length];
        // Only stub undefined methods.
        if (!console[method]) {
            console[method] = noop;
        }
    }

}());

$(document).ready(function() {
    $('input[type=text]#amount_minor_unit').blur(function () {
        $(this).val( ('00' + $(this).val()).substring($(this).val().length, $(this).val().length + 2) );
    });

    // Click event to open the tooltip
    $('a.ringgo-tooltip-click').click(function () {
        // Toggle the use of the show-tooltip class so that we can hide/show it
        $(this).toggleClass('show-tooltip');
    });

    $('a.ringgo-tooltip,a.ringgo-tooltip-click').each(function() {
        // Tooltips on mobile do not always stay on the screen, and therefore require a little nudge in the right direction... No pun intended
        var offsetLeft = 10; // Set the default to 10 so that if the object doesn't exist then it will not crash
        if ($(this).find('.tooltiptext') == undefined) {
            offsetLeft = $(this).find('.tooltiptext').offset().left;
        }
        // Loop to continuously move the tooltip to the right for as long as it is off the screen.
        while(offsetLeft <= 0) {
            // Retrieve the tooltip's left property
            var leftTooltip = $(this).find('.tooltiptext').css('left');
            // As the tooltip is off screen, the offset is negative and therefore nudge it over to the right by 10px
            leftTooltip = parseInt($(this).find('.tooltiptext').css('left')) + 10;
            // Move the tooltip
            $(this).find('.tooltiptext').css('left', leftTooltip + 'px');
            // Retrieve the new offset so that we know when it is now on screen and that we do not loop forever
            offsetLeft = $(this).find('.tooltiptext').offset().left;
        }
    });
    $('ul.error-notification p').each(function() {
        var text = $(this).text();
        var regex = /(.*)(time)(.*?)(passed)(.*?)((\d+)( minutes? ))?((\d+)( seconds?))/gi;
        if(regex.exec(text) !== null) {
            var minutes = parseInt(text.replace(regex, '$7'));
            var seconds = parseInt(text.replace(regex, '$10'));
            var obj = $(this);

            if(!$.isNumeric(minutes)) {
                minutes = 0;
            }

            this.interval = setInterval(function() {
                if(minutes === 0 && seconds === 0) {
                    newText = 'Please try again';
                } else {
                    var newText = '';
                    seconds--;
                    if (seconds < 0) {
                        minutes--;
                        seconds = 59;
                    }
                    if (minutes > 1) {
                        newText = text.replace(regex, '$1$2$3$4$5' + minutes + ' minutes ' + seconds + ' seconds');
                    } else if (minutes === 1) {
                        newText = text.replace(regex, '$1$2$3$4$5' + minutes + ' minute ' + seconds + ' seconds');
                    } else {
                        newText = text.replace(regex, '$1$2$3$4$5' + seconds + ' seconds');
                    }
                }
                $(obj).html(newText);
            }, 1000)
        }
    })
});
// Add Array.indexOf() capability for IE8/9
if (!Array.prototype.indexOf)
{
    Array.prototype.indexOf = function(elt /*, from*/)
    {
        var len = this.length >>> 0;

        var from = Number(arguments[1]) || 0;
        from = (from < 0)
            ? Math.ceil(from)
            : Math.floor(from);
        if (from < 0)
            from += len;

        for (; from < len; from++)
        {
            if (from in this &&
                this[from] === elt)
                return from;
        }
        return -1;
    };
}

function scrollToPos($element, speed){
    $('html, body').animate({ scrollTop: $element.offset().top }, speed);
}


Permits = {
    proofByPost: function (){
        var $proofByPost = $('input[name=PostProof]');

        if($proofByPost.is(':checked')) {
            hideProofs(true);
            $('.proofaddress-wrapper').show(400);
            $('#proofaddress').show(400);
            $('#CurrentProofFieldset').hide(400);
        }else{
            showProofs();
            $('.proofaddress-wrapper').hide(400);
            $('#proofaddress').hide(400);
            $('#CurrentProofFieldset').show(400);
        }

        setTimeout(function(){
            window.scrollTo(0, $proofByPost.offset().top-30);
        },450);
    }
};


window.Utils = {
    
    ajaxMessageRemove: function(name){
        if(name){
            $('.ajaxMessage-'+ name).remove();
        //  Without the param remove all the messages
        } else {
            $('.ajaxMessage').remove();
        }
    },
    
    ajaxMessage: function(message, name, spinner){
        var spinner = spinner || false;
        
        var msg = '';
      
        msg += '<span class="ajaxMessage-'+ name +' ajaxMessage">';
        if(spinner){
            msg += '<img src="/images/ajax-loader.gif" /> ';
        }

        msg += message;
        msg += '</span>';

        return msg;
    }
    
};

/**
 * links acting as pseudo form submits
 * trigger on hyperlink by adding class="pseudo-submit"
 * allow the link to add custom params (eg to signify a different submit action etc) by adding a query string representation
 * to the data-submit="" attribute (eg data-attribute="action=foo&additionalparam=bar")
 */
setupPseudoSubmit = function() {
    $('a.pseudo-submit').on('click', function(e) {
        e.preventDefault();
        var form = $(this).closest('form'),
            params = [];

        if ($(this).data('submit')) {
            params = $.deparam($(this).data('submit'));
            $.each(params, function(name, value) {
                form.append('<input type="hidden" name="' + name + '" value="' + value + '" />');
            });
        }

        form.submit();
    });
}

$(function () {
   setupPseudoSubmit();
});