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
