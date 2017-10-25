
//  For old browsers
if (!window.console) window.console = {log: function() {}};

/*
 * Enable alternative styling when long select options are present
 */
window.selectSizeAdjustment = function selectSizeAdjustment($this){
    var html = new String($this.html());
        if(html.match(/<option([^>]*)>([^<]{60,})<\/option>/gi)){
        $this.parents('li').addClass('singleColumn');
    } else {
        $this.parents('li').removeClass('singleColumn');
    }
};

$(document).ready(function() {
    
    var url = document.URL;
    
    $('#content').delegate('.page-head a.go-up', 'click', function() {
        /*if this is a labyrinth form (i.e. there is a lab id element) use the previous button if there is one*/
        if($('[name="lab_id"]').length && $('.labprev').length) {
            $('.labprev').click();
        } else {
            history.back();
        }
        return false;
    });

    //  Fancy chunky headers
    var $h1 = $('h1')
        .not('#tariffboardtitle h1')
        .not('#tariffboardbgpane h1')
        .not('h1.helpcentre')
    ;
    $h1.wrap('<div class="page-head-wrap"></div>').wrap('<div class="page-head"></div>');
    $h1.before('<a class="go-up icon-arrow-left"><img class="backicon" src="/images/westminster/backarrow.png" alt="RingGo" height="20" width="20" /></a>');

    //  Long dropdowns
    $('.tablelessform select').each(function(){
        selectSizeAdjustment($(this));
    });

    $('span.required').each(function() {
        var innerHtml = $(this).parent().html();
        if(/\s*(\<span class\=\"required\"\>)(\s*\*\s*)(\<\/span\>)(&nbsp;)?/g.test(innerHtml)) {
            innerHtml = innerHtml.replace(/\s*(\<span class\=\"required\"\>)(\s*\*\s*)(\<\/span\>)(&nbsp;)?/g, '\&nbsp\;\<span class="required">\*\<\/span\>');
            $(this).parent().html(innerHtml);
        }
    });
    
    $("#datepicker").datepicker({ altformat: 'd MM yy', 
        changeMonth: true, changeYear: true, 
        yearRange: '-50:+50',
        buttonImage: '/images/datepicker/calendar.gif',
        buttonImageOnly: true,
        showButtonPanel: true,
        showOn: 'both',
        dateFormat: 'd MM yy',
        maxDate: '+1y'});
    $("#datepicker2").datepicker({ altformat: 'd MM yy',
        changeMonth: true, changeYear: true,
        maxDate: '+3m 0y',
        yearRange: '-50:+50',
        buttonImage: '/images/datepicker/calendar.gif', 
        buttonImageOnly: true,
        showButtonPanel: true,
        showOn: 'both',
        dateFormat: 'd MM yy',
        maxDate: '+1y'});
    $("#dateofbirthnopic").datepicker({ altformat: 'dd M yy', 
        changeMonth: true, changeYear: true, 
        yearRange: '-100:+100',
        buttonImage: '/images/datepicker/calendar.gif', 
        buttonImageOnly: true,
        showButtonPanel: false,
        dateFormat: 'dd M yy',
        maxDate: '+0d'});
    $("#dateofbirth").datepicker({ altformat: 'dd M YY',
        changeMonth: true, changeYear: true, 
        yearRange: '-100:+100',
        buttonImage: '/images/datepicker/calendar.gif', 
        buttonImageOnly: true,
        showButtonPanel: true,
        showOn: 'both',
        dateFormat: 'dd M yy',
        maxDate: '+0d'});
    
    


    $("ul.dropdown li").hover(function(){
		
        $(this).addClass("hover");
        $('ul:first',this).css('visibility', 'visible');

    }, function(){

        $(this).removeClass("hover");
        $('ul:first',this).css('visibility', 'hidden');

    });

    $("ul.dropdown li ul li:has(ul)").find("a:first").append(" &raquo; ");



    var hiddenfields = $('.hiddenfieldset');


    var moveAddAddressTextToRight = function() {
        var num = 0;
        // This is to add the class to all of the add address links that should be moving over to the right but haven't been found yet
        $('.wcc-new label > a#addaddress').each(function() {
            $(this).addClass('addAddressLink');
        });
        //currently only runs for WCC. Need to pull out the add address text and float it to the right.
        $('.wcc-new label > .addAddressLink').each(function() {
            num++;

            // Hide the link because for some reason it won't move itself if the replace happens - which needs to happen because BR tags are everywhere
            $(this).hide();
            $(this).addClass('addNew' + num);
            // Get the destination for where the new "fake" link is going to go to
            var newSibing = $(this).closest('li').find("select");

            // This block just removes the br tags because ugh
            var originalParent = $(this).parent();
            var parentHtml = originalParent.html();
            parentHtml = parentHtml.replace(/\s*?\<br.*?\>\s*?/g, '');
            originalParent.html(parentHtml);

            // New "fake" link that will trigger the original through a click event
            var addAddressLink = '<a class="addAddressLink" data-target="addNew' + num + '" id="addAddressLink' + num + '">Add new address</a>';
            //Put the new "fake" link into the destination
            newSibing.after(addAddressLink);

            // This is the triggering click event being added
            $(document).on('click', '#addAddressLink' + num, function() {
                var target = $(this).data('target');
                $('.' + target).click();
            });
        });
    };
    
    
    var moveAddVehicleTextToRight = function() {
        //currently only runs for WCC. Need to pull out the add vehicle text and float it to the right.
        var AddVehicleLinks = $(".wcc-new label > span.addVehicleLink");
        for (var i = 0; i < AddVehicleLinks.length;) {
            var currentLink = $(AddVehicleLinks[i++]);
            var oldParent = currentLink.parent();
            var newParent = $("div.div-vehicleoptions" + i + " > div.elemErrorWrap select");
            newParent.after(currentLink.text("Add registration"));
            var vehiclenumber = oldParent.html().replace(/\D/g,'');
            if (vehiclenumber != "") {
                vehiclenumber = vehiclenumber + ". ";
            }
            oldParent.html('Registration:'+ (i == 1 ? '<span class="required">*</span>' : '') );
            var OwnershipTypeLabel = $("#label-OwnershipId" + i);
            var NewVRNLabel = $("#label-selectvehicle" + i);
            OwnershipTypeLabel.html('Ownership type:'+ (i == 1 ? '<span class="required">*</span>' : '') );
            NewVRNLabel.html('New number plate:'+ (i == 1 ? '<span class="required">*</span>' : '') );
        }
    }
    
    var reformatCheckboxesWithWrappers = function() {
        //This is a hack(ish) so we can style the checkboxes while keeping our origional use of labels intact and work cross-browser.
        
        var updateFakeCheckboxUI = function(elem) {
            //The real checkbox has now been clicked. Lets Update our Fake tick UI to match the state of the real checkbox.
            var tickicon = $($(elem).parent().children()[0]).children()[0];
            if (elem.checked) {
                $(tickicon).attr('style', '');
            } else {
                $(tickicon).attr('style', 'display: none;');
            }
        }
        
        var fakeCheckboxClick = function() {
            //Our fake checkbox has been clicked. Lets find the real checkbox and trigger a click on that too.
            var checkbox = $(this).parent().children()[1];
            $(checkbox).click();
        }
        
        //We need to wrap each checkbox in a container that we can style and observe for clicks on.
        var parents = $('input[type="checkbox"]').parent();
        for (var i = 0; i < parents.length; i++) {
            var parent = $(parents[i]);
            var html = parent.html();
            html = '<div class="checkbox-wrapper"><div class="check-watcher"><div class="checkbox-ticked" style="display: none;"></div></div>' + html + '</div>';
            parent.html(html);
        }
        
        //Wire up the jQ listeners
        $('input[type="checkbox"]').each(function() {
            updateFakeCheckboxUI(this);
        });
        $('input[type="checkbox"]').change(function() {
            updateFakeCheckboxUI(this);
        });
        $('div.check-watcher').on('click', fakeCheckboxClick);
    }
    
    if (url.toLowerCase().indexOf("westminster") >= 0 && $("body.wcc-new").length > 0) {
        //Run these on WCC new
        moveAddAddressTextToRight();
        moveAddVehicleTextToRight();
        reformatCheckboxesWithWrappers();
    }
    
    $('input[name=PostProof]').change(function() {
        proofByPost();
    });

    function proofByPost(){
        var $proofByPost = $('input[name=PostProof]');

        if($proofByPost.is(':checked')) {
            hideProofs();
            $('.proofaddress-wrapper').show(400);
            $('#proofaddress').show(400);
            $('#CurrentProofFieldset').hide(400);
        }else{
            showProofs();
            $('.proofaddress-wrapper').hide(400);
            $('#proofaddress').hide(400);
            $('#CurrentProofFieldset').show(400);
        }
    }

    function hideProofs() {
        speed = 400;
        $('.proofs').closest('li').hide(speed);
        $('.nonTempProofSectionHead').hide(speed);
        $('.postProofDisplay').hide(speed);
    }

    function showProofs() {
        speed = 400;
        $('.proofs').closest('li').show(speed);
        $('.nonTempProofSectionHead').show(speed);
        $('.postProofDisplay').show(speed);
    }
    
});