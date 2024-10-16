document.addEventListener('DOMContentLoaded', function () {
    const backButton = document.querySelector('.back-button');
    backButton.addEventListener('click', function (event) {
        event.preventDefault(); // Prevent default link behavior

        if (document.referrer) {
            // If there is a referrer, go back
            history.back();
        } else {
            // No referrer, redirect to home page
            window.location.href = "{{ url_for('home') }}";
        }
    });
});
