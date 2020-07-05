// Javascript: Toggles radiobutton checked setting.
function toggle(eID, boolSetting) {
    eID.checked = boolSetting;
}

// Javascript: Reveals elements that were hidden (by opacity = 0) before vcf load.
function revealElements(classN) {
    els = document.getElementsByClassName(classN);
    for (var i = 0, n = els.length; i < n; i++) {
        document.getElementsByClassName(classN)[i].style.opacity = 1;
    }
}

// Javascript: Hides elements that by setting opacity = 0.
function hideElements(classN) {
    els = document.getElementsByClassName(classN);
    for (var i = 0, n = els.length; i < n; i++) {
        document.getElementsByClassName(classN)[i].style.opacity = 0;
    }
}

// Javascript: Bog standard checkbox group toggle.
function toggleGroup(eIDstring, boolSetting) {
    checkboxes = document.getElementsByName(eIDstring);
    for (var i = 0, n = checkboxes.length; i < n; i++) {
        toggle(checkboxes[i], boolSetting);
    }
}


// Javascript: Standard checkbox group toggle with shiny twist.
function shinyToggleGroup(eIDstring, boolSetting) {
    toggleGroup(eIDstring, boolSetting)

    // Now set up a string to pass back to Shiny to confirm the action.
    if (boolSetting)
        rtnString = 'SET_ALL';
    else
        rtnString = 'CLEAR_ALL';

    // Make sure it propogates up to the R shiny event handler.
    // I was finding that some events were missed without this line..
    Shiny.setInputValue(eIDstring, rtnString, {
        priority: 'event'
    });
}


// Message handler to change the backbround image of an element when triggered from Rshiny.
// (used to add loading notification etc.)
Shiny.addCustomMessageHandler("changeBackgroundImage_h_id", changeBackgroundImage);

function changeBackgroundImage(args) {
    document.getElementById(args.eName).style.backgroundImage = "url(" + args.img + ")"
}


// Message handler to hide UI elements while load in progress.
Shiny.addCustomMessageHandler("hideUIelementsOnLoad_h_id", hideUIelementsOnLoad);

function hideUIelementsOnLoad(className) {
    hideElements(className);
}


// Message handler to reveal UI elements once everything's loaded.
Shiny.addCustomMessageHandler("revealUIelementsOnLoad_h_id", revealUIelementsOnLoad);

function revealUIelementsOnLoad(className) {
    revealElements(className);
}

